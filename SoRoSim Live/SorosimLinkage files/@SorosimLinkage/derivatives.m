%Function to calculate qdot and qdotdot at a given time used for odex function of matlab
%Last modified by Anup Teejo Mathew 13.06.2022
function dqdt=derivatives(Tr,t,qqd,uqt)
%lsqoptions is used if custom actuator strength is enabled and actuator strength is calculated using lsqlin or lsqnonneg
% 
persistent tlast
if t==0
    tlast=cputime;
end
if cputime-tlast>0.5
 		tlast = cputime;
        disp(t);
end

ndof   = Tr.ndof;
q      = qqd(1:ndof);
qd     = qqd(ndof+1:2*ndof);

n_jact = Tr.n_jact;
N      = Tr.N;
nsig   = Tr.nsig;

if ~Tr.CAS
    if Tr.Actuated
        WrenchControlled = Tr.WrenchControlled;
        i_jactq          = Tr.i_jactq;
        for i=1:n_jact
            if ~WrenchControlled(i)
                q(i_jactq(i))  = uqt{i}{1}(t);
                qd(i_jactq(i)) = uqt{i}{2}(t);
            end
        end
    end
end

M = zeros(ndof,ndof);
C = zeros(ndof,ndof);
F = zeros(ndof,1);

g   = zeros(4*nsig,4);
J   = zeros(6*nsig,ndof);
eta = zeros(6*nsig,1);
Jd  = zeros(6*nsig,ndof); % for dependent base Jd is Z

G = Tr.G; %gravity

dof_start = 1; %starting dof of current piece
i_sig     = 1;
i_sig_nj  = 1;
g_ini     = Tr.g_ini; %initial configuration of all link wrt its previous link
g_Ltip    = repmat(eye(4),N,1);
J_Ltip    = repmat(zeros(6,ndof),N,1);
Jd_Ltip   = repmat(zeros(6,ndof),N,1); 

eta_Ltip   = zeros(N*6,1); %total velocity J*qd+eta_t

iLpre     = Tr.iLpre;

for i=1:N

    if iLpre(i)>0
        g_here       = g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        J_here       = Ad_g_ini_inv*J_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        Jd_here      = Ad_g_ini_inv*Jd_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        
        eta_here   = Ad_g_ini_inv*eta_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6);
    else
        g_here   = g_ini((i-1)*4+1:i*4,:);
        J_here   = zeros(6,ndof);
        Jd_here  = zeros(6,ndof);
        
        eta_here   = zeros(6,1);
    end
    
    %Joint
    dof_here = Tr.CVTwists{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    qd_here  = qd(dof_start:dof_start+dof_here-1);
    B_here   = Tr.CVTwists{i}(1).B;
    xi_star  = Tr.CVTwists{i}(1).xi_star;

    if dof_here==0 %fixed joint (N)
        g_joint   = eye(4);
        TgB_here  = zeros(6,ndof);
        TgBd_here = zeros(6,ndof);
    else
        xi          = B_here*q_here+xi_star;
        xid         = B_here*qd_here;

        [g_joint,Tg,Tgd]=variable_expmap_gTgTgd_mex(xi,xid);

        TgB_here                                    = zeros(6,ndof);
        TgB_here(:,dof_start:dof_start+dof_here-1)  = Tg*B_here;
        TgBd_here                                   = zeros(6,ndof);
        TgBd_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*Tg*B_here+Tgd*B_here;
    end

    %updating g, Jacobian, Jacobian_dot and eta
    g_here         = g_here*g_joint;
    Ad_g_joint_inv = dinamico_Adjoint(ginv(g_joint));
    J_here         = Ad_g_joint_inv*(J_here+TgB_here);
    Jd_here        = Ad_g_joint_inv*(Jd_here+TgBd_here);
    eta_here       = Ad_g_joint_inv*(eta_here+TgB_here(:,dof_start:dof_start+dof_here-1)*qd_here);

    g((i_sig-1)*4+1:i_sig*4,:)  = g_here;
    J((i_sig-1)*6+1:i_sig*6,:)  = J_here;
    Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
    eta((i_sig-1)*6+1:i_sig*6)  = eta_here;
    i_sig                       = i_sig+1;

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'

        gi        = Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here    = g_here*gi;
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        J_here    = Ad_gi_inv*J_here;
        Jd_here   = Ad_gi_inv*Jd_here;
        eta_here  = Ad_gi_inv*eta_here;

        g((i_sig-1)*4+1:i_sig*4,:)  = g_here;
        J((i_sig-1)*6+1:i_sig*6,:)  = J_here;
        Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
        eta((i_sig-1)*6+1:i_sig*6)  = eta_here;
        
        M_here = Tr.VLinks(Tr.LinkIndex(i)).M;
        if Tr.Gravity
            F = F+J_here'*M_here*dinamico_Adjoint(ginv(g_here))*G;
        end
        
        if ~isempty(Tr.M_added)
            M_here_add = M_here+Tr.M_added((i_sig_nj-1)*6+1:i_sig_nj*6,:); 
        else
            M_here_add = M_here;
        end

        Qtemp  = J_here'*M_here_add; % temporary variable to avoid repetetion
        
        M     = M+Qtemp*J_here;
        C     = C+Qtemp*Jd_here+J_here'*dinamico_coadj(eta_here)*M_here*J_here;
        
        % bringing all quantities to the end of rigid link
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here    = g_here*gf;
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        J_here    = Ad_gf_inv*J_here;
        Jd_here   = Ad_gf_inv*Jd_here;
        eta_here  = Ad_gf_inv*eta_here;
        
        i_sig    = i_sig+1;
        i_sig_nj = i_sig_nj+1;
    end

    dof_start = dof_start+dof_here;

    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link

        dof_here   = Tr.CVTwists{i}(j+1).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        qd_here    = qd(dof_start:dof_start+dof_here-1);
        
        gi      = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        ld      = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        
        if ~Tr.CVTwists{i}(j+1).CI %not custom integration
            xi_star = Tr.CVTwists{i}(j+1).xi_star;
            Ms      = Tr.CVTwists{i}(j+1).Ms;
            Xs      = Tr.CVTwists{i}(j+1).Xs;
            Ws      = Tr.CVTwists{i}(j+1).Ws;
            nip     = Tr.CVTwists{i}(j+1).nip;
        else
            [nip,Xs,Ws] = Tr.CVTwists{i}(j+1).CIFn(Tr,i,j,t,q_here);
            Tr.CVTwists{i}(j+1).nip = nip;
            Tr.CVTwists{i}(j+1).Xs = Xs;
            Tr.CVTwists{i}(j+1).Ws = Ws;
            Tr.CVTwists{i}(j+1).Updatexi_star();
            Tr.CVTwists{i}(j+1).UpdateMEG();
            xi_star = Tr.CVTwists{i}(j+1).xi_star;
            Ms      = Tr.CVTwists{i}(j+1).Ms;
        end

        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here    = g_here*gi;
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        J_here    = Ad_gi_inv*J_here;
        Jd_here   = Ad_gi_inv*Jd_here;
        eta_here  = Ad_gi_inv*eta_here;

        %scaling of quantities using the formula: Lscale m = 1 unit
        
        Lscale         = ld;
        g_here(1:3,4)  = g_here(1:3,4)/Lscale;
        J_here(4:6,:)  = J_here(4:6,:)/Lscale;
        Jd_here(4:6,:) = Jd_here(4:6,:)/Lscale;
        eta_here(4:6)  = eta_here(4:6)/Lscale;
        G = G/Lscale;

        g((i_sig-1)*4+1:i_sig*4,:)  = g_here;
        J((i_sig-1)*6+1:i_sig*6,:)  = J_here;
        Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
        eta((i_sig-1)*6+1:i_sig*6)  = eta_here;
       
        ii = 1;
        if Ws(ii)>0
            W_here  = Ws(ii);
            Ms_here = Ms(6*(ii-1)+1:6*ii,:);

            if isempty(Tr.M_added)
                Ms_here_add = Ms_here;
            else
                Ms_here_add = Ms_here+Tr.M_added((i_sig_nj-1)*6+1:i_sig_nj*6,:);
            end

            %scaling
            Ms_here(1:3,:) = Ms_here(1:3,:)/Lscale;
            Ms_here(4:6,:) = Ms_here(4:6,:)*Lscale;
            if Tr.Gravity
                F = F+W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %rescale
            end

            Ms_here_add(1:3,:) = Ms_here_add(1:3,:)/Lscale;
            Ms_here_add(4:6,:) = Ms_here_add(4:6,:)*Lscale;

            Qtemp = W_here*J_here'*Ms_here_add*Lscale^2;%rescale (temporary variable to avoid repeated computation)
            M     = M+Qtemp*J_here; 
            C     = C+Qtemp*Jd_here+W_here*J_here'*dinamico_coadj(eta_here)*Ms_here*J_here*Lscale^2; %rescale
        end
        
        i_sig    = i_sig+1;
        i_sig_nj = i_sig_nj+1;

        for ii=2:nip

            H = Xs(ii)-Xs(ii-1);
            
            if Tr.Z_order==4
                
                xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1),2); 
                xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
                xi_Z1here(1:3) = xi_Z1here(1:3)*Lscale; %scaling
                xi_Z2here(1:3) = xi_Z2here(1:3)*Lscale;
                
                %B is Phi for independent basis
                    
                B_Z1here  = Tr.CVTwists{i}(j+1).B_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
                B_Z2here  = Tr.CVTwists{i}(j+1).B_Z2(6*(ii-2)+1:6*(ii-1),:);

                if dof_here>0

                    xi_Z1here = B_Z1here*q_here+xi_Z1here;
                    xi_Z2here = B_Z2here*q_here+xi_Z2here;

                    xid_Z1here  = B_Z1here*qd_here;
                    ad_xi_Z1here = dinamico_adj(xi_Z1here);

                    BGamma_here  = (H/2)*(B_Z1here+B_Z2here)+... %dBqdq = B
                                       ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);

                    Gammadd_Z4_dq_here = ((sqrt(3)*H^2)/6)*dinamico_adj(xid_Z1here)*B_Z2here; 
                    Gammad_here   = BGamma_here*qd_here;

                else

                    ad_xi_Z1here = dinamico_adj(xi_Z1here);
                    BGamma_here  = (H/2)*(B_Z1here+B_Z2here)+...
                                       ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);
                    Gammadd_Z4_dq_here = zeros(6,dof_here); 
                    Gammad_here   = zeros(6,1); 

                end
                Gamma_here = (H/2)*(xi_Z1here+xi_Z2here)+...
                             ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;
                      
            else % order 2
                
                xi_Zhere = xi_star(6*(ii-2)+1:6*(ii-1),4);
                xi_Zhere(1:3) = xi_Zhere(1:3)*Lscale; %scaling
                    
                B_Zhere  = Tr.CVTwists{i}(j+1).B_Z(6*(ii-2)+1:6*(ii-1),:);%note this step

                if dof_here>0
                    xi_Zhere = B_Zhere*q_here+xi_Zhere;
                    BGamma_here = H*B_Zhere;
                    Gammad_here     = BGamma_here*qd_here;
                else
                    BGamma_here = H*B_Zhere;
                    Gammad_here   = zeros(6,1); 
                end
                Gamma_here  = H*xi_Zhere;

            end

            [gh,TGamma_here,TGammad_here] = variable_expmap_gTgTgd_mex(Gamma_here,Gammad_here); % mex code, C program

            TBGamma_here                                    = zeros(6,ndof);
            TBGamma_here(:,dof_start:dof_start+dof_here-1)  = TGamma_here*BGamma_here;
            TBGammad_here                                   = zeros(6,ndof);
            TBGammad_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*TBGamma_here(:,dof_start:dof_start+dof_here-1)+TGammad_here*BGamma_here;

            if Tr.Z_order==4
                TBGammad_here(:,dof_start:dof_start+dof_here-1) = TBGammad_here(:,dof_start:dof_start+dof_here-1)+TGamma_here*Gammadd_Z4_dq_here;
            end
            
            %updating g, Jacobian, Jacobian_dot and eta
            g_here     = g_here*gh;
            Ad_gh_inv  = dinamico_Adjoint(ginv(gh));
            J_here     = Ad_gh_inv*(J_here+TBGamma_here); %full
            Jd_here    = Ad_gh_inv*(Jd_here+TBGammad_here); %full
            eta_here   = Ad_gh_inv*(eta_here+TBGamma_here(:,dof_start:dof_start+dof_here-1)*qd_here);

            
            g((i_sig-1)*4+1:i_sig*4,:)  = g_here;
            J((i_sig-1)*6+1:i_sig*6,:)  = J_here;
            Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
            eta((i_sig-1)*6+1:i_sig*6)  = eta_here;
            
            %integrals evaluation
            if Ws(ii)>0
                W_here  = Ws(ii);
                Ms_here = Ms(6*(ii-1)+1:6*ii,:);
                
                if isempty(Tr.M_added)
                    Ms_here_add = Ms_here;
                else
                    Ms_here_add = Ms_here+Tr.M_added((i_sig_nj-1)*6+1:i_sig_nj*6,:);
                end
                 
                %scaling
                Ms_here(1:3,:) = Ms_here(1:3,:)/Lscale;
                Ms_here(4:6,:) = Ms_here(4:6,:)*Lscale;
                if Tr.Gravity
                    F = F+W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %rescale
                end
                
                Ms_here_add(1:3,:) = Ms_here_add(1:3,:)/Lscale;
                Ms_here_add(4:6,:) = Ms_here_add(4:6,:)*Lscale;
                
                Qtemp = W_here*J_here'*Ms_here_add*Lscale^2;%rescale (temporary variable to avoid repeated computation)
                M     = M+Qtemp*J_here; 
                C     = C+Qtemp*Jd_here+W_here*J_here'*dinamico_coadj(eta_here)*Ms_here*J_here*Lscale^2; %rescale
            end
            i_sig    = i_sig+1;
            i_sig_nj = i_sig_nj+1;
        end

        %scaling back quantities
        g_here(1:3,4)  = g_here(1:3,4)*Lscale;
        J_here(4:6,:)  = J_here(4:6,:)*Lscale;
        Jd_here(4:6,:) = Jd_here(4:6,:)*Lscale;
        eta_here(4:6)  = eta_here(4:6)*Lscale;
        
        G              = G*Lscale;

        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
        g_here    = g_here*gf;
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        J_here    = Ad_gf_inv*J_here;
        Jd_here   = Ad_gf_inv*Jd_here;
        eta_here  = Ad_gf_inv*eta_here;

        dof_start = dof_start+dof_here;
    end
    
    g_Ltip((i-1)*4+1:i*4,:)    = g_here;
    J_Ltip((i-1)*6+1:i*6,:)    = J_here;
    Jd_Ltip((i-1)*6+1:i*6,:)   = Jd_here;
    eta_Ltip((i-1)*6+1:i*6,:)  = eta_here;

end

%% Point Force
if Tr.PointForce
    F = F+ComputePointForce(Tr,J,g,t);
end

%% Joint and soft link actuation
if Tr.Actuated

    nact = Tr.nact;
    Bq   = zeros(ndof,nact);
    u    = zeros(nact,1);

    %revolute, prismatic, helical joints
    n_jact          = Tr.n_jact;
    Bqj1            = Tr.Bqj1;
    n_1dof          = size(Bqj1,2);
    Bq(:,1:n_1dof)  = Tr.Bqj1;

    %for other joints
    i_jact          = Tr.i_jact;
    i_u             = n_1dof+1;
    i_jactq         = Tr.i_jactq;

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

    %Cable actuation
    n_sact=Tr.n_sact;

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

    if ~Tr.CAS % if Custom Actuation Strength u is found in the next section
        for i=1:n_jact 
            if ~WrenchControlled(i) %swapping for joint cooridnate controlled joints, u becomes qdd_joint
                u(i)             = uqt{i}{3}(t);
            else
                u(i)             = uqt{i}(t);
            end
        end
        for i=n_jact+1:nact
            u(i) = uqt{i}(t);
        end
    end
else
    Bq=0;
    u=0;
end

%% Custom External Force, Actuator or Actuator Strength
if Tr.CEFP||Tr.CAP||Tr.CAS

    i_sig    = 1;
    g_act    = g;
    J_act    = J;
    Jd_act   = Jd;
    eta_act  = eta;
    
    for i=1:N
        i_sig = i_sig+1;
        if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
            i_sig = i_sig+1;
        end
        for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
            Lscale = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
            for ii=1:Tr.CVTwists{i}(j+1).nip
                g_here        = g((i_sig-1)*4+1:i_sig*4,:);
                g_here(1:3,4) = g_here(1:3,4)*Lscale;
                g_act((i_sig-1)*4+1:i_sig*4,:) = g_here;

                J_here        = J((i_sig-1)*6+1:i_sig*6,:);
                J_here(4:6,:) = J_here(4:6,:)*Lscale;
                J_act((i_sig-1)*6+1:i_sig*6,:) = J_here;
                
                eta_here      = eta((i_sig-1)*6+1:i_sig*6);
                eta_here(4:6) = eta_here(4:6,:)*Lscale;
                eta_act((i_sig-1)*6+1:i_sig*6) = eta_here;

                Jd_here        = Jd((i_sig-1)*6+1:i_sig*6,:);
                Jd_here(4:6,:) = Jd_here(4:6,:)*Lscale;
                Jd_act((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
                i_sig = i_sig+1;               
            end
        end
    end
    
    M_act   = M;
    C_act   = C;
    F_act   = F;
    Bq_act  = Bq;
    
    if Tr.CEFP
        Fext = CustomExtForce(Tr,q,g_act,J_act,t,qd,eta_act,Jd_act);
        FextP = CustomExtPointForce(Tr,q,g_act,J_act,t,qd,eta_act,Jd_act);
        F_custom = CustomExtForce_Qspace(Tr,J,Fext,FextP);
        F = F+F_custom;
    end
    if Tr.CAP
        Fact = CustomActuation(Tr,q,g_act,J_act,t,qd,eta_act,Jd_act);
        TauC = CustomActuation_Qspace(Tr,Fact,q,t);
        F    = F+TauC; %added with F
    end
    if Tr.CAS  %customized actuator strength (eg. Controller)
        u = CustomActuatorStrength(Tr,q,g_act,J_act,t,qd,eta_act,Jd_act,M_act,C_act,F_act,Bq_act);
    end

end

%% Closed Loop Joints
if Tr.nCLj>0

    A  = zeros(Tr.CLprecompute.nCLp,Tr.ndof);
    Ad = zeros(Tr.CLprecompute.nCLp,Tr.ndof);
    e  = zeros(Tr.CLprecompute.nCLp,1);

    k=1;
    for ii=1:Tr.nCLj

        Bp     = Tr.CLprecompute.Bp{ii};
        i_sigA = Tr.CLprecompute.i_sigA(ii);
        i_sigB = Tr.CLprecompute.i_sigB(ii);

        if Tr.iACL(ii)>0
            LinkA = Tr.VLinks(Tr.LinkIndex(Tr.iACL(ii)));
            gA    = g((i_sigA-1)*4+1:i_sigA*4,:);
            JA    = J((i_sigA-1)*6+1:i_sigA*6,:);
            JdA   = Jd((i_sigA-1)*6+1:i_sigA*6,:);

            if LinkA.linktype=='s'
                gA(1:3,4)  = gA(1:3,4)*LinkA.lp{end};
                JA(4:6,:)  = JA(4:6,:)*LinkA.lp{end};
                JdA(4:6,:) = JdA(4:6,:)*LinkA.lp{end};
                
                gf  = LinkA.gf{end};
            else
                gf  = LinkA.gf;
            end
            
            gA  = gA*gf;
            Ad_gf_inv=dinamico_Adjoint(ginv(gf));
            JA  = Ad_gf_inv*JA;
            JdA = Ad_gf_inv*JdA;
        else
            gA  = eye(4);
            JA  = zeros(6,Tr.ndof);
            JdA = zeros(6,Tr.ndof);
        end

        if Tr.iCLB(ii)>0
            LinkB = Tr.VLinks(Tr.LinkIndex(Tr.iCLB(ii)));
            gB    = g((i_sigB-1)*4+1:i_sigB*4,:);
            JB    = J((i_sigB-1)*6+1:i_sigB*6,:);
            JdB   = Jd((i_sigB-1)*6+1:i_sigB*6,:);

            if LinkB.linktype=='s'
                gB(1:3,4)  = gB(1:3,4)*LinkB.lp{end};
                JB(4:6,:)  = JB(4:6,:)*LinkB.lp{end};
                JdB(4:6,:) = JdB(4:6,:)*LinkB.lp{end};
                
                gf  = LinkB.gf{end};
            else
                gf  = LinkB.gf;
            end
            gB  = gB*gf;
            Ad_gf_inv = dinamico_Adjoint(ginv(gf));
            JB  = Ad_gf_inv*JB;
            JdB = Ad_gf_inv*JdB;
        else
            gB  = eye(4);
            JB  = zeros(6,Tr.ndof);
            JdB = zeros(6,Tr.ndof);
        end
        
        gCLjA = gA*Tr.gACLj{ii};
        gCLjB = gB*Tr.gBCLj{ii};
        JA    = dinamico_Adjoint(ginv(Tr.gACLj{ii}))*JA;
        JB    = dinamico_Adjoint(ginv(Tr.gBCLj{ii}))*JB; %moving to CLj frame

        gCLjAB = ginv(gCLjA)*gCLjB;
        Ad_gCLjAB_inv = dinamico_Adjoint(ginv(gCLjAB));
        JA     = Ad_gCLjAB_inv*JA;
        JdA    = Ad_gCLjAB_inv*JdA;

        A(k:k+size(Bp,2)-1,:)  = Bp'*(JA-JB); %change
        Ad(k:k+size(Bp,2)-1,:) = Bp'*(JdA-JdB); %change
        e(k:k+size(Bp,2)-1,:)  = Bp'*piecewise_logmap(ginv(gCLjAB));

        k=k+size(Bp,2);
    end

    e(e.*e<1e-12)=0;
    if rank(A)<size(A,1)
        [A,iRows] = LIRows(A);
        Ad        = Ad(iRows,:);
        e         = e(iRows);
    end

    if e'*e<1e-12
        e = zeros(size(e));
    end
end
%%

if Tr.Damped
    D = Tr.D;
else
    D = 0;
end
K = Tr.K;

if Tr.nCLj>0
    T   = Tr.T_BS;

    if rank(M)<Tr.ndof %M singular
        nA = size(A,1);
        MA  = [M A';A zeros(nA,nA)];
        CDA = [C+D;Ad+(2/T)*A];
        KA  = [K;zeros(nA,Tr.ndof)];
        BqA  = [Bq;zeros(nA,Tr.nact)];
        FA  = [F;-(1/T^2)*e];

        for i=1:n_jact
            if ~WrenchControlled(i) %swapping for joint cooridnate controlled joints
                MA_temp           = MA; %shift after concatnation
                MA(:,i_jactq(i))  = -BqA(:,i); 
                BqA(:,i)          = -MA_temp(:,i_jactq(i));  %Tau=Bq*u
            end
        end

        if rank(MA)<(Tr.ndof+Tr.CLprecompute.nCLp) %MA singular
            qdd_lambda = pinv(MA)*(BqA*u+FA-KA*q-CDA*qd);
        else
            qdd_lambda = MA\(BqA*u+FA-KA*q-CDA*qd);
        end
        qdd = qdd_lambda(1:Tr.ndof); %lambda is discarded
    else
        for i=1:n_jact 
            if ~WrenchControlled(i) %swapping for joint cooridnate controlled joints
                M_temp           = M;
                M(:,i_jactq(i))  = -Bq(:,i); 
                Bq(:,i)          = -M_temp(:,i_jactq(i));  %Tau=Bq*u
            end
        end
        P   = eye(Tr.ndof)-A'*(A*M^-1*A')^-1*A*M^-1;
        qdd = M\(P*(Bq*u+F-K*q-(C+D)*qd)-A'*(A*M^-1*A')^-1*(Ad*qd+(2/T)*A*qd+(1/T^2)*e));
    end
else
    for i=1:n_jact
        if ~WrenchControlled(i)  %swapping for joint cooridnate controlled joints
            M_temp           = M;
            M(:,i_jactq(i))  = -Bq(:,i); 
            Bq(:,i)          = -M_temp(:,i_jactq(i));  %Tau=Bq*u
        end
    end
    if rank(M)<Tr.ndof
        qdd = pinv(M)*(Bq*u+F-K*q-(C+D)*qd);
    else
        qdd = M\(Bq*u+F-K*q-(C+D)*qd);
    end
end



if Tr.Actuated
    for i=1:n_jact
        if ~WrenchControlled(i)
            qdd(i_jactq(i)) = u(i); %replace u with qdd if wrench controlled
        end
    end
end

dqdt=[qd;qdd]; %implicit time
end