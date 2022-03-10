%Function to calculate qdot and qdotdot at a given time used for ode45 dynamics
%Last modified by Anup Teejo Mathew 02.03.2022
function dqdt=derivatives(Tr,t,qqd,uqt)

persistent tlast
if t==0
    tlast=cputime;
end
if cputime-tlast>0.5
 		tlast = cputime;
        disp(t);
end

M_scale = Tr.q_scale*Tr.q_scale';
%precomputed values
if Tr.Damped
    D = Tr.D.*M_scale;
else
    D = 0;
end
K = Tr.K.*M_scale;

ndof   = Tr.ndof;
q      = qqd(1:ndof);
qd     = qqd(ndof+1:2*ndof);
n_jact = Tr.n_jact;
N      = Tr.N;
nsig   = Tr.nsig;

if Tr.Actuated
    WrenchControlled = Tr.WrenchControlled;
    i_jactq          = Tr.i_jactq;
    for i=1:n_jact
        if ~WrenchControlled(i)
            q(i_jactq(i))   = uqt{i}{1}(t);
            qd(i_jactq(i))  = uqt{i}{2}(t);
        end
    end
end

M = zeros(ndof,ndof);
C = zeros(ndof,ndof);
F = zeros(ndof,1);
G = Tr.G;

g   = zeros(4*nsig,4);
J   = zeros(6*nsig,ndof);
Jd  = zeros(6*nsig,ndof);
eta = zeros(6*nsig,1);

dof_start = 1; %starting dof of current piece
i_sig     = 1;
g_ini     = Tr.g_ini;
g_Ltip    = repmat(eye(4),N,1);
J_Ltip    = repmat(zeros(6,ndof),N,1);
Jd_Ltip   = repmat(zeros(6,ndof),N,1);
iLpre     = Tr.iLpre;

for i=1:N

    if iLpre(i)>0
        g_here   = g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        J_here   = J_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        J_here   = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)))*J_here;
        Jd_here  = Jd_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        Jd_here  = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)))*Jd_here;
        eta_here = J_here*qd;
    else
        g_here   = g_ini((i-1)*4+1:i*4,:);
        J_here   = zeros(6,ndof);
        Jd_here  = zeros(6,ndof);
        eta_here = J_here*qd;
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
        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='U' %special case for universal joint. Considered as 2 revolute joints
            % first revolute joint
            xi          = B_here(:,1)*q_here(1)+xi_star;
            xid         = B_here(:,1)*qd_here(1);
            theta_here  = norm(xi(1:3));
            thetad_here = (xid(1:3)'*xi(1:3))/theta_here;
            g_joint     = joint_expmap(xi);

            Tg                    = variable_Texpmap(1,theta_here,xi);
            TgB_here              = zeros(6,ndof);
            TgB_here(:,dof_start) = Tg*B_here(:,1);

            Tgd                    = variable_dotTexpmap(1,theta_here,thetad_here,xi,xid);
            TgBd_here              = zeros(6,ndof);
            TgBd_here(:,dof_start) = dinamico_adj(eta_here)*Tg*B_here(:,1)+Tgd*B_here(:,1);

            g_here     = g_here*g_joint;
            J_here     = dinamico_Adjoint(ginv(g_joint))*(J_here+TgB_here);
            Jd_here    = dinamico_Adjoint(ginv(g_joint))*(Jd_here+TgBd_here);
            eta_here   = J_here*qd;

            % second revolute joint
            xi          = B_here(:,2)*q_here(2)+xi_star;
            xid         = B_here(:,2)*qd_here(2);
            theta_here  = norm(xi(1:3));
            thetad_here = (xid(1:3)'*xi(1:3))/theta_here;
            g_joint     = joint_expmap(xi);

            Tg                      = variable_Texpmap(1,theta_here,xi);
            TgB_here                = zeros(6,ndof);
            TgB_here(:,dof_start+1) = Tg*B_here(:,2);

            Tgd                      = variable_dotTexpmap(1,theta_here,thetad_here,xi,xid);
            TgBd_here                = zeros(6,ndof);
            TgBd_here(:,dof_start+1) = dinamico_adj(eta_here)*Tg*B_here(:,2)+Tgd*B_here(:,2);
        else
            xi          = B_here*q_here+xi_star;
            xid         = B_here*qd_here;
            theta_here  = norm(xi(1:3));
            thetad_here = (xid(1:3)'*xi(1:3))/theta_here;
            g_joint     = joint_expmap(xi);

            Tg                                         = variable_Texpmap(1,theta_here,xi);
            TgB_here                                   = zeros(6,ndof);
            TgB_here(:,dof_start:dof_start+dof_here-1) = Tg*B_here;

            Tgd                                         = variable_dotTexpmap(1,theta_here,thetad_here,xi,xid);
            TgBd_here                                   = zeros(6,ndof);
            TgBd_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*Tg*B_here+Tgd*B_here;
        end
    end

    %updating g, Jacobian, Jacobian_dot and eta
    g_here     = g_here*g_joint;
    J_here     = dinamico_Adjoint(ginv(g_joint))*(J_here+TgB_here);
    Jd_here    = dinamico_Adjoint(ginv(g_joint))*(Jd_here+TgBd_here);
    eta_here   = J_here*qd;

    g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
    J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
    Jd((i_sig-1)*6+1:i_sig*6,:)   = Jd_here;
    eta((i_sig-1)*6+1:i_sig*6)    = eta_here;
    i_sig                         = i_sig+1;

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'

        gi       =  Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here   = g_here*gi;
        J_here   = dinamico_Adjoint(ginv(gi))*J_here;
        Jd_here  = dinamico_Adjoint(ginv(gi))*Jd_here;
        eta_here = J_here*qd;

        g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
        J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
        Jd((i_sig-1)*6+1:i_sig*6,:)   = Jd_here;
        eta((i_sig-1)*6+1:i_sig*6)    = eta_here;
        i_sig                         = i_sig+1;

        M_here = Tr.VLinks(Tr.LinkIndex(i)).Ms;
        M      = M+J_here'*M_here*J_here;
        C      = C+J_here'*(M_here*Jd_here+dinamico_coadj(eta_here)*M_here*J_here);
        if Tr.Gravity
            F  = F+J_here'*M_here*dinamico_Adjoint(ginv(g_here))*G;
        end

        % bringing all quantities to the end of rigid link
        gf       =  Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here   = g_here*gf;
        J_here   = dinamico_Adjoint(ginv(gf))*J_here;
        Jd_here  = dinamico_Adjoint(ginv(gf))*Jd_here;
        eta_here = J_here*qd;
    end

    dof_start = dof_start+dof_here;

    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link

        dof_here   = Tr.CVTwists{i}(j+1).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        qd_here    = qd(dof_start:dof_start+dof_here-1);
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
        Jd_here       = dinamico_Adjoint(ginv(gi))*Jd_here;
        eta_here      = J_here*qd;

        %scaling of quantities
        Lscale         = ld;
        g_here(1:3,4)  = g_here(1:3,4)/Lscale;
        J_here(4:6,:)  = J_here(4:6,:)/Lscale;
        Jd_here(4:6,:) = Jd_here(4:6,:)/Lscale;
        eta_here(4:6)  = eta_here(4:6)/Lscale;
        ld             = 1;
        G              = G/Lscale;


        g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
        J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
        Jd((i_sig-1)*6+1:i_sig*6,:)   = Jd_here;
        eta((i_sig-1)*6+1:i_sig*6)    = eta_here;
        i_sig                         = i_sig+1;

        for ii=2:nGauss

            B_Z1here      = B_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
            B_Z2here      = B_Z2(6*(ii-2)+1:6*(ii-1),:);
            xi_starZ1here = xi_star(6*(ii-2)+1:6*(ii-1),2);
            xi_starZ2here = xi_star(6*(ii-2)+1:6*(ii-1),3);

            xi_starZ1here(1:3) = xi_starZ1here(1:3)*Lscale; %scaling
            xi_starZ2here(1:3) = xi_starZ2here(1:3)*Lscale;

            xi_Z1here     = B_Z1here*q_here+xi_starZ1here;
            xi_Z2here     = B_Z2here*q_here+xi_starZ2here;

            H             = (Xs(ii)-Xs(ii-1))*ld;
            Gamma_here    = (H/2)*(xi_Z1here+xi_Z2here)+...
                ((sqrt(3)*H^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
            k_here        = Gamma_here(1:3);
            theta_here    = norm(k_here);
            gh            = variable_expmap(theta_here,Gamma_here);

            BGamma_here   = (H/2)*(B_Z1here+B_Z2here)+...
                ((sqrt(3)*H^2)/12)*(dinamico_adj(xi_Z1here)*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);

            Gammad_here   = BGamma_here*qd_here;
            kd_here       = Gammad_here(1:3);
            thetad_here   = (kd_here'*k_here)/theta_here;

            TGamma_here                                    = variable_Texpmap(1,theta_here,Gamma_here);
            TBGamma_here                                   = zeros(6,ndof);
            TBGamma_here(:,dof_start:dof_start+dof_here-1) = TGamma_here*BGamma_here;



            TGammad_here                                    = variable_dotTexpmap(1,theta_here,thetad_here,Gamma_here,Gammad_here);
            TBGammad_here                                   = zeros(6,ndof);
            TBGammad_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*TGamma_here*BGamma_here+TGammad_here*BGamma_here;

            %updating g, Jacobian, Jacobian_dot and eta
            g_here     = g_here*gh;
            J_here     = dinamico_Adjoint(ginv(gh))*(J_here+TBGamma_here);
            Jd_here    = dinamico_Adjoint(ginv(gh))*(Jd_here+TBGammad_here);
            eta_here   = J_here*qd;

            g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
            J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
            Jd((i_sig-1)*6+1:i_sig*6,:)   = Jd_here;
            eta((i_sig-1)*6+1:i_sig*6)    = eta_here;
            i_sig                         = i_sig+1;

            %integrals evaluation
            if ii<nGauss
                W_here   = Ws(ii);
                Ms_here  = Ms(6*(ii-1)+1:6*ii,:);
                %scaling
                Ms_here(1:3,:)   = Ms_here(1:3,:)/Lscale;
                Ms_here(4:6,:)   = Ms_here(4:6,:)*Lscale;
                M        = M+ld*W_here*J_here'*Ms_here*J_here*Lscale^2; %rescale
                C        = C+ld*W_here*J_here'*(Ms_here*Jd_here+dinamico_coadj(eta_here)*Ms_here*J_here)*Lscale^2; %rescale
                if Tr.Gravity
                    F    = F+ld*W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %rescale
                end
            end

        end

        %scaling back quantities
        g_here(1:3,4)  = g_here(1:3,4)*Lscale;
        J_here(4:6,:)  = J_here(4:6,:)*Lscale;
        Jd_here(4:6,:) = Jd_here(4:6,:)*Lscale;
        G              = G*Lscale;

        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gf       = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
        g_here   = g_here*gf;
        J_here   = dinamico_Adjoint(ginv(gf))*J_here;
        Jd_here  = dinamico_Adjoint(ginv(gf))*Jd_here;
        eta_here = J_here*qd;

        dof_start = dof_start+dof_here;
    end
    g_Ltip((i-1)*4+1:i*4,:)  = g_here;
    J_Ltip((i-1)*6+1:i*6,:)  = J_here;
    Jd_Ltip((i-1)*6+1:i*6,:) = Jd_here;
end

%% Point Force
if Tr.PointForce
    F = F+ComputePointForce(Tr,J,t,g);
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


    n_ljact=length(i_jact);
    for iii=i_u:n_ljact

        i=i_jact(iii);

        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='U'
            Bq(i_jactq(i_u:i_u+1),i_u:i_u+1) = [1 0;0 1];
            i_u = i_u+2;
        elseif Tr.VLinks(Tr.LinkIndex(i)).jointtype=='C'
            Bq(i_jactq(i_u:i_u+1),i_u:i_u+1) = [1 0;0 1];
            i_u = i_u+2;
        else

            i_sig = 1;
            for ii=1:i-1
                i_sig = i_sig+1;
                for jj=1:Tr.VLinks(Tr.LinkIndex(ii)).npie-1
                    i_sig = i_sig+1+Tr.VLinks(Tr.LinkIndex(ii)).nGauss{jj};
                end
                if Tr.VLinks(Tr.LinkIndex(ii)).linktype=='r'
                    i_sig = i_sig+1;
                end
            end

            J_here = J((i_sig-1)*6+1:i_sig*6,:);
            S_here = J_here(:,i_jactq(i_u:i_u+Tr.CVTwists{i}(1).dof-1));
            B_here = Tr.CVTwists{i}(1).B;

            Bq(i_jactq(i_u:i_u+Tr.CVTwists{i}(1).dof-1),i_u:i_u+Tr.CVTwists{i}(1).dof-1) = S_here'*B_here;
            i_u = i_u+Tr.CVTwists{i}(1).dof;

        end
    end

    %Cable actuation
    n_sact=Tr.n_sact;

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
                M_temp           = M;
                M(:,i_jactq(i))  = Bq(:,i);
                Bq(:,i)          = M_temp(:,i_jactq(i));  %Tau=Bq*u
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

    J_scale  = repmat(Tr.q_scale',6,1);
    i_sig    = 1;
    g_act    = zeros(nsig*4,4);
    J_act    = zeros(nsig*6,ndof);
    eta_act  = zeros(nsig*6,1);
    Jd_act   = zeros(nsig*6,ndof);
    for i=1:N
        % joint
        g_act((i_sig-1)*4+1:i_sig*4,:)  = g((i_sig-1)*4+1:i_sig*4,:);
        J_act((i_sig-1)*6+1:i_sig*6,:)  = J((i_sig-1)*6+1:i_sig*6,:)./J_scale;
        eta_act((i_sig-1)*6+1:i_sig*6)  = eta((i_sig-1)*6+1:i_sig*6);
        Jd_act((i_sig-1)*6+1:i_sig*6,:) = Jd((i_sig-1)*6+1:i_sig*6,:)./J_scale;
        i_sig = i_sig+1;
        if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
            g_act((i_sig-1)*4+1:i_sig*4,:)  = g((i_sig-1)*4+1:i_sig*4,:);
            J_act((i_sig-1)*6+1:i_sig*6,:)  = J((i_sig-1)*6+1:i_sig*6,:)./J_scale;
            eta_act((i_sig-1)*6+1:i_sig*6)  = eta((i_sig-1)*6+1:i_sig*6);
            Jd_act((i_sig-1)*6+1:i_sig*6,:) = Jd((i_sig-1)*6+1:i_sig*6,:)./J_scale;
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

                eta_here      = eta((i_sig-1)*6+1:i_sig*6);
                eta_here(4:6) = eta_here(4:6)*Lscale;
                eta_act((i_sig-1)*6+1:i_sig*6)  = eta_here;

                Jd_here        = Jd((i_sig-1)*6+1:i_sig*6,:);
                Jd_here(4:6,:) = Jd_here(4:6,:)*Lscale;
                Jd_act((i_sig-1)*6+1:i_sig*6,:)  = Jd_here./J_scale;
                i_sig = i_sig+1;
            end
        end
    end
    M_act   = M./M_scale;
    C_act   = C./M_scale;
    F_act   = F./Tr.q_scale;
    B_scale = repmat(Tr.q_scale,1,Tr.nact);
    Bq_act  = Bq./B_scale;
    
    if Tr.CEFP
        Fext = CustomExtForce(Tr,q.*Tr.q_scale,g_act,J_act,t,qd.*Tr.q_scale,eta_act,Jd_act);
        [F_custom,M_custom,C_custom] = CustomExtForce_Qspace(Tr,J,Fext,Jd,eta);
        M = M+M_custom;
        F = F+F_custom;
        C = C+C_custom;
    end
    if Tr.CAP
        Fact = CustomActuation(Tr,q.*Tr.q_scale,g_act,J_act,t,qd.*Tr.q_scale,eta_act,Jd_act);
        TauC = CustomActuation_Qspace(Tr,Fact);
        F    = F+TauC; %added with F
    end
    if Tr.CAS
        u = CustomActuatorStrength(Tr,q.*Tr.q_scale,g_act,J_act,t,qd.*Tr.q_scale,eta_act,Jd_act,M_act,C_act,F_act,Bq_act);
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
                gA  = gA*gf;
                JA  = dinamico_Adjoint(ginv(gf))*JA;
                JdA = dinamico_Adjoint(ginv(gf))*JdA;
            else
                gf  = LinkA.gf;
                gA  = gA*gf;
                JA  = dinamico_Adjoint(ginv(gf))*JA;
                JdA = dinamico_Adjoint(ginv(gf))*JdA;
            end
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
                gB  = gB*gf;
                JB  = dinamico_Adjoint(ginv(gf))*JB;
                JdB = dinamico_Adjoint(ginv(gf))*JdB;
            else
                gf  = LinkB.gf;
                gB  = gB*gf;
                JB  = dinamico_Adjoint(ginv(gf))*JB;
                JdB = dinamico_Adjoint(ginv(gf))*JdB;
            end
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
        JA     = dinamico_Adjoint(ginv(gCLjAB))*JA;
        JdA    = dinamico_Adjoint(ginv(gCLjAB))*JdA;

        A(k:k+size(Bp,2)-1,:)  = Bp'*(JA-JB); %change
        Ad(k:k+size(Bp,2)-1,:) = Bp'*(JdA-JdB); %change
        e(k:k+size(Bp,2)-1,:)  = Bp'*piecewise_logmap(1,ginv(gCLjAB));

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

if Tr.nCLj>0
    P   = eye(Tr.ndof)-A'*(A*M^-1*A')^-1*A*M^-1;
    T   = Tr.T_BS;
    qdd = M\(P*(Bq*u+F-K*q-(C+D)*qd)-A'*(A*M^-1*A')^-1*(Ad*qd+(2/T)*A*qd+(1/T^2)*e));
else
    qdd = M\(Bq*u+F-K*q-(C+D)*qd);
end


if Tr.Actuated
    for i=1:n_jact
        if ~WrenchControlled(i)
            qdd(i_jactq(i)) = u(i);
        end
    end
end

dqdt=[qd;qdd];
end