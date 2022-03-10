%Evaluates the static equilibrium equation for fsolve
%Last modified by Anup Teejo Mathew 02.03.2022
function E=Equilibrium(Tr,qu,uq,magnifier)

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
K_scale    = Tr.q_scale*Tr.q_scale';
K          = Tr.K.*K_scale;

J          = zeros(6*n_sig,ndof);
g          = zeros(4*n_sig,4);
dof_start  = 1; %starting dof of current piece
i_sig      = 1;
F          = zeros(ndof,1);%external force
g_ini      = Tr.g_ini;
g_Ltip     = repmat(eye(4),N,1);
J_Ltip     = repmat(zeros(6,ndof),N,1);
iLpre      = Tr.iLpre;
M 	       = zeros(ndof,ndof);


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
        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='U' %special case for universal joint. Considered as 2 revolute joints
            % first revolute joint
            xi         = B_here(:,1)*q_here(1)+xi_star;
            theta_here = norm(xi(1:3));
            g_joint    = joint_expmap(xi);

            Tg                    = variable_Texpmap(1,theta_here,xi);
            TgB_here              = zeros(6,ndof);
            TgB_here(:,dof_start) = Tg*B_here(:,1);

            g_here     = g_here*g_joint;
            J_here     = dinamico_Adjoint(ginv(g_joint))*(J_here+TgB_here);

            % second revolute joint
            xi         = B_here(:,2)*q_here(2)+xi_star;
            theta_here = norm(xi(1:3));
            g_joint    = joint_expmap(xi);

            Tg                      = variable_Texpmap(1,theta_here,xi);
            TgB_here                = zeros(6,ndof);
            TgB_here(:,dof_start+1) = Tg*B_here(:,2);
        else
            xi         = B_here*q_here+xi_star;
            theta_here = norm(xi(1:3));
            g_joint    = joint_expmap(xi);

            Tg                                         = variable_Texpmap(1,theta_here,xi);
            TgB_here                                   = zeros(6,ndof);
            TgB_here(:,dof_start:dof_start+dof_here-1) = Tg*B_here;
        end
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

        M_here = Tr.VLinks(Tr.LinkIndex(i)).Ms;
        M      = M+J_here'*M_here*J_here;
        if Tr.Gravity
            F = F+J_here'*M_here*dinamico_Adjoint(ginv(g_here))*G;
        end

        % bringing all quantities to the end of rigid link
        gf     = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here = g_here*gf;
        J_here = dinamico_Adjoint(ginv(gi))*J_here;
    end

    dof_start = dof_start+dof_here;

    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link

        dof_here   = Tr.CVTwists{i}(j+1).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        B_Z1       = Tr.CVTwists{i}(j+1).B_Z1;
        B_Z2       = Tr.CVTwists{i}(j+1).B_Z2;
        xi_star    = Tr.CVTwists{i}(j+1).xi_star;
        gi         = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        ld         = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Ms         = Tr.VLinks(Tr.LinkIndex(i)).Ms{j};
        Xs         = Tr.VLinks(Tr.LinkIndex(i)).Xs{j};
        Ws         = Tr.VLinks(Tr.LinkIndex(i)).Ws{j};
        nGauss     = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};

        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here        = g_here*gi;
        J_here        = dinamico_Adjoint(ginv(gi))*J_here;

        %scaling of quantities
        Lscale        = ld;
        g_here(1:3,4) = g_here(1:3,4)/Lscale;
        J_here(4:6,:) = J_here(4:6,:)/Lscale;
        ld            = 1;
        G             = G/Lscale;

        g((i_sig-1)*4+1:i_sig*4,:) = g_here;
        J((i_sig-1)*6+1:i_sig*6,:) = J_here;
        i_sig                      = i_sig+1;

        for ii=2:nGauss

            B_Z1here      = B_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
            B_Z2here      = B_Z2(6*(ii-2)+1:6*(ii-1),:);

            xi_starZ1here = xi_star(6*(ii-2)+1:6*(ii-1),2);
            xi_starZ2here = xi_star(6*(ii-2)+1:6*(ii-1),3);

            xi_starZ1here(1:3) = xi_starZ1here(1:3)*Lscale; %scaling
            xi_starZ2here(1:3) = xi_starZ2here(1:3)*Lscale;

            xi_Z1here     = B_Z1here*q_here+xi_starZ1here;
            xi_Z2here     = B_Z2here*q_here+xi_starZ2here;

            h             = (Xs(ii)-Xs(ii-1))*ld;
            Gamma_here    = (h/2)*(xi_Z1here+xi_Z2here)+...
                            ((sqrt(3)*h^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
            k_here        = Gamma_here(1:3);
            theta_here    = norm(k_here);
            gh            = variable_expmap(theta_here,Gamma_here);

            BGamma_here   = (h/2)*(B_Z1here+B_Z2here)+...
                            ((sqrt(3)*h^2)/12)*(dinamico_adj(xi_Z1here)*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);

            TGamma_here                                    = variable_Texpmap(1,theta_here,Gamma_here);
            TBGamma_here                                   = zeros(6,ndof);
            TBGamma_here(:,dof_start:dof_start+dof_here-1) = TGamma_here*BGamma_here;

            %updating g, Jacobian, Jacobian_dot and eta
            g_here     = g_here*gh;
            J_here     = dinamico_Adjoint(ginv(gh))*(J_here+TBGamma_here);

            g((i_sig-1)*4+1:i_sig*4,:)  = g_here;
            J((i_sig-1)*6+1:i_sig*6,:)  = J_here;
            i_sig                       = i_sig+1;

            %integrals evaluation
            if ii<nGauss
                W_here  = Ws(ii);
                Ms_here = Ms(6*(ii-1)+1:6*ii,:);
                %scaling
                Ms_here(1:3,:)   = Ms_here(1:3,:)/Lscale;
                Ms_here(4:6,:)   = Ms_here(4:6,:)*Lscale;
                M                = M+ld*W_here*J_here'*Ms_here*J_here*Lscale^2; %rescale
                if Tr.Gravity
                    F = F+(ld)*W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %scaled back for addition
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
    F = F+ComputePointForce(Tr,J,0,g);
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
    i_jact           = Tr.i_jact;
    i_tau            = n_1dof+1;
    i_jactq          = Tr.i_jactq;

    for i=i_jact(i_tau:end)
        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='U'
            Bq(i_jactq(i_tau:i_tau+1),i_tau:i_tau+1) = [1 0;0 1];
            i_tau = i_tau+2;
        elseif Tr.VLinks(Tr.LinkIndex(i)).jointtype=='C'
            Bq(i_jactq(i_tau:i_tau+1),i_tau:i_tau+1) = [1 0;0 1];
            i_tau = i_tau+2;
        elseif Tr.VLinks(Tr.LinkIndex(i)).jointtype=='A'
            i_sig = 1;
            f     = 1;
            for ii=1:i-1
                i_sig = i_sig+1;
                for jj=1:Tr.VLinks(ii).npie-1
                    i_sig = i_sig+1+Tr.VLinks(ii).nGauss{jj};
                end
                if Tr.VLinks(ii).linktype=='r'
                    i_sig = i_sig+1;
                end
                f=f+Tr.VLinks(ii).npie;
            end
            J_here = J((i_sig-1)*6+1:i_sig*6,:);
            S_here = J_here(:,i_jactq(i_tau:i_tau+2));
            B_here = Tr.CVTwists{f}.B;

            Bq(i_jactq(i_tau:i_tau+2),i_tau:i_tau+2) = S_here'*B_here;
            i_tau = i_tau+3;
        elseif Tr.VLinks(Tr.LinkIndex(i)).jointtype=='S'
            i_sig = 1;
            for ii=1:i-1
                i_sig = i_sig+1;
                for jj=1:Tr.VLinks(ii).npie-1
                    i_sig = i_sig+1+Tr.VLinks(ii).nGauss{jj};
                end
                if Tr.VLinks(ii).linktype=='r'
                    i_sig = i_sig+1;
                end
            end
            J_here = J((i_sig-1)*6+1:i_sig*6,:);
            S_here = J_here(:,i_jactq(i_tau:i_tau+2));
            B_here = [eye(3);zeros(3,3)];

            Bq(i_jactq(i_tau:i_tau+2),i_tau:i_tau+2) = S_here'*B_here;
            i_tau = i_tau+3;
        else %free joint
            i_sig = 1;
            for ii=1:i-1
                i_sig = i_sig+1;
                for jj=1:Tr.VLinks(ii).npie-1
                    i_sig = i_sig+1+Tr.VLinks(ii).nGauss{jj};
                end
                if Tr.VLinks(ii).linktype=='r'
                    i_sig = i_sig+1;
                end
            end
            J_here = J((i_sig-1)*6+1:i_sig*6,:);
            S_here = J_here(:,i_jactq(i_tau:i_tau+5));
            B_here = eye(6);

            Bq(i_jactq(i_tau:i_tau+5),i_tau:i_tau+5) = S_here'*B_here;
            i_tau = i_tau+6;
        end
    end

    %cable actuation
    n_sact = Tr.n_sact;

    for ii=1:n_sact

        dcii = cell(1,N); dcpii = cell(1,N); Sdivii = cell(1,N); Edivii = cell(1,N);

        for i=1:N
            dcii{i}   = Tr.dc{ii,i};
            dcpii{i}  = Tr.dcp{ii,i};
            Sdivii{i} = Tr.Sdiv{ii,i};
            Edivii{i} = Tr.Ediv{ii,i};
        end
        Insideii      = Tr.Inside{ii};

        if Insideii
            Bq(:,n_jact+ii) = ComputeCableActuation(Tr,dcii,dcpii,Sdivii,Edivii,q);
        else
            Bq(:,n_jact+ii) = ComputeCableActuation2(Tr,dcii,Sdivii,Edivii,J,g);
        end

    end
    if ~Tr.CAS
        WrenchControlled = Tr.WrenchControlled;
        for i=1:n_jact
            if ~WrenchControlled(i)
                u(i) = qu(i_jactq(i));
            else
                u(i) = uq(i);
            end
        end

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

    J_scale  = repmat(Tr.q_scale',6,1);
    i_sig    = 1;
    g_act    = zeros(n_sig*4,4);
    J_act    = zeros(n_sig*6,ndof);

    for i=1:N
        % joint
        g_act((i_sig-1)*4+1:i_sig*4,:)  = g((i_sig-1)*4+1:i_sig*4,:);
        J_act((i_sig-1)*6+1:i_sig*6,:)  = J((i_sig-1)*6+1:i_sig*6,:)./J_scale;
        i_sig = i_sig+1;
        if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
            g_act((i_sig-1)*4+1:i_sig*4,:)  = g((i_sig-1)*4+1:i_sig*4,:);
            J_act((i_sig-1)*6+1:i_sig*6,:)  = J((i_sig-1)*6+1:i_sig*6,:)./J_scale;
            i_sig = i_sig+1;
        end
        for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
            Lscale = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
            for ii=1:Tr.VLinks(Tr.LinkIndex(i)).nGauss{j}

                g_here        = g((i_sig-1)*4+1:i_sig*4,:);
                g_here(1:3,4) = g_here(1:3,4)*Lscale;
                g_act((i_sig-1)*4+1:i_sig*4,:)  = g_here;

                J_here        = J((i_sig-1)*6+1:i_sig*6,:);
                J_here(4:6,:) = J_here(4:6,:)*Lscale;
                J_act((i_sig-1)*6+1:i_sig*6,:)  = J_here./J_scale;

                i_sig = i_sig+1;
            end
        end
    end

    if Tr.CEFP
        Fext = CustomExtForce(Tr,q.*Tr.q_scale,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof));
        [F_custom,M_custom,~] = CustomExtForce_Qspace(Tr,J,Fext,zeros(6*n_sig,ndof),zeros(6*n_sig,1));
        M = M+M_custom;
        F = F+F_custom;
    end
    if Tr.CAP
        Fact = CustomActuation(Tr,q.*Tr.q_scale,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof));
        TauC = CustomActuation_Qspace(Tr,Fact);
        F    = F+TauC; %added with F
    end
    if Tr.CAS
        u = CustomActuatorStrength(Tr,q.*Tr.q_scale,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof));
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
        e(k:k+size(Bp,2)-1,:) = Bp'*piecewise_logmap(1,ginv(gCLjAB));

        k = k+size(Bp,2);
    end

    if rank(A)<size(A,1)
        [A,iRows] = LIRows(A);
        e         = e(iRows); %iRows are the linearly independent rows of A
    end

end

%% Equilibrium
if Tr.nCLj>0
    P       = eye(Tr.ndof)-A'*(A*M^-1*A')^-1*A*M^-1;
    E_Force = LIRows(P)*(Bq*u+F-K*q);
    E       = [E_Force;e];
else
    E = K*q-Bq*u-F;
end

if magnifier
    E=E*1e7;
end

end