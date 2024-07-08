%Evaluates the static equilibrium equation for fsolve
%Last modified by Anup Teejo Mathew 02.03.2022
function E=Equilibrium(Tr,qu,uq,magnifier,lsqoptions)

N      = Tr.N;
ndof   = Tr.ndof;
n_sig  = Tr.nsig;
G      = Tr.G;
q      = qu;

if Tr.Actuated
    n_jact           = Tr.n_jact;
    WrenchControlled = Tr.WrenchControlled;
    i_jactq          = Tr.i_jactq;
    for i=1:n_jact
        if ~WrenchControlled(i)
            q(i_jactq(i))   = uq(i);
        end
    end
end


J          = zeros(6*n_sig,ndof);
g          = zeros(4*n_sig,4);
dof_start  = 1; %starting dof of current piece
i_sig      = 1;
F          = zeros(ndof,1);%external force
g_ini      = Tr.g_ini;
g_Ltip     = repmat(eye(4),N,1);
J_Ltip     = repmat(zeros(6,ndof),N,1);
iLpre      = Tr.iLpre;
% M 	       = zeros(ndof,ndof);


for i=1:N

    if iLpre(i)>0
        g_here = g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        J_here = J_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        J_here = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)))*J_here;
    else
        g_here    = g_ini((i-1)*4+1:i*4,:);
        J_here    = zeros(6,ndof);
    end

    %Joint
    dof_here   = Tr.CVTwists{i}(1).dof;
    q_here     = q(dof_start:dof_start+dof_here-1);
    B_here     = Tr.CVTwists{i}(1).B;
    xi_star    = Tr.CVTwists{i}(1).xi_star;

    if dof_here==0 %fixed joint (N)
        g_joint     = eye(4);
        TgB_here    = zeros(6,ndof);
    else
        
        xi           = B_here*q_here+xi_star;
        [g_joint,Tg] = variable_expmap_gTg(xi); %mex file is slightly slower

        TgB_here                                   = zeros(6,ndof);
        TgB_here(:,dof_start:dof_start+dof_here-1) = Tg*B_here;

    end

    %updating g, Jacobian, Jacobian_dot and eta
    g_here     = g_here*g_joint;
    J_here     = dinamico_Adjoint(ginv(g_joint))*(J_here+TgB_here);

    g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
    J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
    i_sig                         = i_sig+1;

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'

        gi     = Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here = g_here*gi;
        J_here = dinamico_Adjoint(ginv(gi))*J_here;

        g((i_sig-1)*4+1:i_sig*4,:) = g_here;
        J((i_sig-1)*6+1:i_sig*6,:) = J_here;
        i_sig                      = i_sig+1;

        M_here = Tr.VLinks(Tr.LinkIndex(i)).M;
        if Tr.Gravity
            F = F+J_here'*M_here*dinamico_Adjoint(ginv(g_here))*G;
        end

        % bringing all quantities to the end of rigid link
        gf     = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here = g_here*gf;
        J_here = dinamico_Adjoint(ginv(gf))*J_here;
    end

    dof_start = dof_start+dof_here;

    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link

        dof_here   = Tr.CVTwists{i}(j+1).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        
        gi      = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        ld      = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        
        if ~Tr.CVTwists{i}(j+1).CI %not custom integration
            xi_star = Tr.CVTwists{i}(j+1).xi_star;
            Ms      = Tr.CVTwists{i}(j+1).Ms;
            Xs      = Tr.CVTwists{i}(j+1).Xs;
            Ws      = Tr.CVTwists{i}(j+1).Ws;
            nip     = Tr.CVTwists{i}(j+1).nip;
        else
            [nip,Xs,Ws] = Tr.CVTwists{i}(j+1).CIFn(Tr,i,j,0,q_here);
            Tr.CVTwists{i}(j+1).nip = nip;
            Tr.CVTwists{i}(j+1).Xs = Xs;
            Tr.CVTwists{i}(j+1).Ws = Ws;
            Tr.CVTwists{i}(j+1).Updatexi_star();
            Tr.CVTwists{i}(j+1).UpdateMEG();
            xi_star = Tr.CVTwists{i}(j+1).xi_star;
            Ms      = Tr.CVTwists{i}(j+1).Ms;
        end


        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here        = g_here*gi;
        J_here        = dinamico_Adjoint(ginv(gi))*J_here;

        %scaling of quantities
        Lscale        = ld;
        g_here(1:3,4) = g_here(1:3,4)/Lscale;
        J_here(4:6,:) = J_here(4:6,:)/Lscale;
        G             = G/Lscale;

        g((i_sig-1)*4+1:i_sig*4,:) = g_here;
        J((i_sig-1)*6+1:i_sig*6,:) = J_here;
        
        ii = 1;
        if Ws(ii)>0
            W_here  = Ws(ii);
            Ms_here = Ms(6*(ii-1)+1:6*ii,:);
            %scaling
            Ms_here(1:3,:)   = Ms_here(1:3,:)/Lscale;
            Ms_here(4:6,:)   = Ms_here(4:6,:)*Lscale;
            if Tr.Gravity
                F = F+W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %scaled back for addition
            end
        end
        
        
        i_sig  = i_sig+1;
        
        for ii=2:nip
            
            H = Xs(ii)-Xs(ii-1);
            
            if Tr.Z_order==4
                
                xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1),2);
                xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
                xi_Z1here(1:3) = xi_Z1here(1:3)*Lscale; %scaling
                xi_Z2here(1:3) = xi_Z2here(1:3)*Lscale;
                
                B_Z1here  = Tr.CVTwists{i}(j+1).B_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
                B_Z2here  = Tr.CVTwists{i}(j+1).B_Z2(6*(ii-2)+1:6*(ii-1),:);
                if dof_here>0
                    xi_Z1here = B_Z1here*q_here+xi_Z1here;
                    xi_Z2here = B_Z2here*q_here+xi_Z2here;
                end
                ad_xi_Z1here = dinamico_adj(xi_Z1here);
                BGamma_here  = (H/2)*(B_Z1here+B_Z2here)+...
                           ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here); 

                Gamma_here   = (H/2)*(xi_Z1here+xi_Z2here)+...
                                ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;
                      
            else % order 2
                xi_Zhere = xi_star(6*(ii-2)+1:6*(ii-1),4);
                xi_Zhere(1:3) = xi_Zhere(1:3)*Lscale; %scaling


                B_Zhere  = Tr.CVTwists{i}(j+1).B_Z(6*(ii-2)+1:6*(ii-1),:);%note this step
                if dof_here>0
                    xi_Zhere = B_Zhere*q_here+xi_Zhere;
                end
                BGamma_here = H*B_Zhere;

                
                Gamma_here  = H*xi_Zhere;

            end
                        
            [gh,TGamma_here] = variable_expmap_gTg(Gamma_here); %mex file is slightly slower

            TBGamma_here                                   = zeros(6,ndof);
            TBGamma_here(:,dof_start:dof_start+dof_here-1) = TGamma_here*BGamma_here;

            %updating g, Jacobian, Jacobian_dot and eta
            g_here     = g_here*gh;
            J_here     = dinamico_Adjoint(ginv(gh))*(J_here+TBGamma_here);

            g((i_sig-1)*4+1:i_sig*4,:)  = g_here;
            J((i_sig-1)*6+1:i_sig*6,:)  = J_here;
            i_sig                       = i_sig+1;

            %integrals evaluation
            if Ws(ii)>0
                W_here  = Ws(ii);
                Ms_here = Ms(6*(ii-1)+1:6*ii,:);
                %scaling
                Ms_here(1:3,:)   = Ms_here(1:3,:)/Lscale;
                Ms_here(4:6,:)   = Ms_here(4:6,:)*Lscale;
                if Tr.Gravity
                    F = F+W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %scaled back for addition
                end
            end

        end

        %scaling back quantities
        g_here(1:3,4) = g_here(1:3,4)*Lscale;
        J_here(4:6,:) = J_here(4:6,:)*Lscale;
        G             = G*Lscale;

        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gf            = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
        g_here        = g_here*gf;
        J_here        = dinamico_Adjoint(ginv(gf))*J_here;

        dof_start = dof_start+dof_here;
    end
    g_Ltip((i-1)*4+1:i*4,:) = g_here;
    J_Ltip((i-1)*6+1:i*6,:) = J_here;
end
%% Point Force
if Tr.PointForce
    F = F+ComputePointForce(Tr,J,g,0);
end

%% Actuation

if Tr.Actuated

    nact = Tr.nact;
    Bq   = zeros(ndof,nact);
    u    = zeros(nact,1);

    %revolute, prismatic, helical joints
    n_jact           = Tr.n_jact;
    Bqj1             = Tr.Bqj1;
    n_1dof           = size(Bqj1,2);
    Bq(:,1:n_1dof)   = Tr.Bqj1;

    %for other joints
    i_jact  = Tr.i_jact;
    i_u     = n_1dof+1;
    i_jactq = Tr.i_jactq;
    
    n_ljact = length(i_jact);
    for iii=n_1dof+1:n_ljact
        i = i_jact(iii);
        i_sig = 1;
        dof_start = 1;
        for ii=1:i-1
            i_sig = i_sig+1;
            dof_start = dof_start+Tr.CVTwists{ii}(1).dof; %joint
            for j=1:Tr.VLinks(Tr.LinkIndex(ii)).npie-1
                i_sig = i_sig+Tr.CVTwists{ii}(j+1).nip;
                dof_start = dof_start+Tr.CVTwists{ii}(j+1).dof;
            end
            if Tr.VLinks(Tr.LinkIndex(ii)).linktype=='r'
                i_sig = i_sig+1;
            end
        end

        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='C'
            dof_here = dof_start+Tr.CVTwists{i}(1).dof;
            Bq(i_jactq(i_u:i_u+dof_here-1),i_u:i_u+dof_here-1) = [1 0;0 1];
            i_u = i_u+2;
        else 
            dof_here = Tr.CVTwists{i}(1).dof;
            J_here = J((i_sig-1)*6+1:i_sig*6,:);
            S_here = J_here(:,i_jactq(i_u:i_u+dof_here-1));
            B_here = Tr.CVTwists{i}(1).B;

            Bq(i_jactq(i_u:i_u+dof_here-1),i_u:i_u+dof_here-1) = S_here'*B_here;
            i_u = i_u+dof_here;
        end
    end


    %cable actuation
    n_sact = Tr.n_sact;

    for ii=1:n_sact

        dcii = cell(1,N); dcpii = cell(1,N); Sdivii = zeros(N,1); Edivii = zeros(N,1);

        for i=1:N
            dcii{i}   = Tr.dc{ii,i};
            dcpii{i}  = Tr.dcp{ii,i};
            Sdivii(i) = Tr.Sdiv(ii,i);
            Edivii(i) = Tr.Ediv(ii,i);
        end
        Insideii      = Tr.Inside(ii);

        if Insideii
            Bq(:,n_jact+ii) = ComputeCableActuation(Tr,dcii,dcpii,Sdivii,Edivii,q);
        else
            Bq(:,n_jact+ii) = ComputeCableActuation2(Tr,dcii,Sdivii,Edivii,J,g);
        end

    end
    
    for i=1:n_jact
        if ~WrenchControlled(i)
            u(i) = qu(i_jactq(i));
        else
            u(i) = uq(i);
        end
    end
    
    if ~Tr.CAS
        for i=n_jact+1:nact
            u(i) = uq(i);
        end
    end
else
    Bq=0;
    u=0;
end
%% Custom External Force, Actuation or Actuator Strength
if Tr.CEFP||Tr.CAP||Tr.CAS

%     J_scale  = repmat(Tr.q_scale',6,1);
    i_sig    = 1;
    g_act    = g;
    J_act    = J;

    for i=1:N
        % joint
        i_sig = i_sig+1;
        if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
            i_sig = i_sig+1;
        end
        for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
            Lscale = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
            for ii=1:Tr.CVTwists{i}(j+1).nip

                g_here        = g((i_sig-1)*4+1:i_sig*4,:);
                g_here(1:3,4) = g_here(1:3,4)*Lscale;
                g_act((i_sig-1)*4+1:i_sig*4,:)  = g_here;

                J_here        = J((i_sig-1)*6+1:i_sig*6,:);
                J_here(4:6,:) = J_here(4:6,:)*Lscale;
                J_act((i_sig-1)*6+1:i_sig*6,:)  = J_here;%./J_scale;

                i_sig = i_sig+1;
            end
        end
    end
    

    F_act   = F;%
    Bq_act  = Bq;%

    if Tr.CEFP
        Fext = CustomExtForce(Tr,q,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof));
        FextP = CustomExtPointForce(Tr,q,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof));
        F_custom = CustomExtForce_Qspace(Tr,J,Fext,FextP);
        F = F+F_custom;
    end
    if Tr.CAP
        Fact = CustomActuation(Tr,q,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof));
        TauC = CustomActuation_Qspace(Tr,Fact,q,0);
        F    = F+TauC; %added with F
    end
    if Tr.CAS
        u = CustomActuatorStrength(Tr,q,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof),zeros(ndof,ndof),zeros(ndof,ndof),F_act,Bq_act,lsqoptions);
    end
end


%% Closed Loop Joints
if Tr.nCLj>0

    A = zeros(Tr.CLprecompute.nCLp,Tr.ndof);
    e = zeros(Tr.CLprecompute.nCLp,1);

    k=1;
    for ii=1:Tr.nCLj

        Bp     = Tr.CLprecompute.Bp{ii};
        i_sigA = Tr.CLprecompute.i_sigA(ii);
        i_sigB = Tr.CLprecompute.i_sigB(ii);

        if Tr.iACL(ii)>0
            LinkA = Tr.VLinks(Tr.LinkIndex(Tr.iACL(ii)));
            gA    = g((i_sigA-1)*4+1:i_sigA*4,:);
            JA    = J((i_sigA-1)*6+1:i_sigA*6,:);

            if LinkA.linktype=='s'
                gA(1:3,4) = gA(1:3,4)*LinkA.lp{end}; %unscale
                JA(4:6,:) = JA(4:6,:)*LinkA.lp{end};
                gf = LinkA.gf{end};
                gA = gA*gf;
                JA = dinamico_Adjoint(ginv(gf))*JA;
            else
                gf = LinkA.gf;
                gA = gA*gf;
                JA = dinamico_Adjoint(ginv(gf))*JA;
            end
        else
            gA = eye(4);
            JA = zeros(6,Tr.ndof);
        end

        if Tr.iCLB(ii)>0
            LinkB = Tr.VLinks(Tr.LinkIndex(Tr.iCLB(ii)));
            gB    = g((i_sigB-1)*4+1:i_sigB*4,:);
            JB    = J((i_sigB-1)*6+1:i_sigB*6,:);

            if LinkB.linktype=='s'
                gB(1:3,4) = gB(1:3,4)*LinkB.lp{end};
                JB(4:6,:) = JB(4:6,:)*LinkB.lp{end};
                gf        = LinkB.gf{end};
                gB        = gB*gf;
                JB        = dinamico_Adjoint(ginv(gf))*JB;
            else
                gf = LinkB.gf;
                gB = gB*gf;
                JB = dinamico_Adjoint(ginv(gf))*JB;
            end
        else
            gB = eye(4);
            JB = zeros(6,Tr.ndof);
        end

        gCLjA = gA*Tr.gACLj{ii};
        gCLjB = gB*Tr.gBCLj{ii};
        JA    = dinamico_Adjoint(ginv(Tr.gACLj{ii}))*JA;
        JB    = dinamico_Adjoint(ginv(Tr.gBCLj{ii}))*JB; %moving to CLj frame

        gCLjAB = ginv(gCLjA)*gCLjB;
        JA     = dinamico_Adjoint(ginv(gCLjAB))*JA; %Transforming Jacobian to appropriate frame

        A(k:k+size(Bp,2)-1,:) = Bp'*(JA-JB);
        e(k:k+size(Bp,2)-1,:) = Bp'*piecewise_logmap(ginv(gCLjAB));

        k = k+size(Bp,2);
    end

    if rank(A)<size(A,1)
        [A,iRows] = LIRows(A);
        e         = e(iRows); %iRows are the linearly independent rows of A. Comment this for some
    end
    if e'*e<1e-12
        e = zeros(size(e));
    end

end

%% Equilibrium

K  = Tr.K;
if Tr.nCLj>0
    P = eye(Tr.ndof)-A'*(A*A')^-1*A;

    E = LIRows(P)*(Bq*u+F-K*q); % P projects away A'*Lambda
%         E = P*(Bq*u+F-K*q);% comment above and uncomment this for some
    E = [E;e];
else
    E = K*q-Bq*u-F;
end


if magnifier
    E=E*1e7;
end

end