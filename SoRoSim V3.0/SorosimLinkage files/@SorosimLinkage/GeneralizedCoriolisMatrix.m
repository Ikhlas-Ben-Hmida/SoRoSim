%Function that calculates the generalized coriolis matrix  at every significant points (24.05.2021)

function C = GeneralizedCoriolisMatrix(Tr,q,qd)

ndof      = Tr.ndof;
C         = zeros(ndof,ndof);
N         = Tr.N;

if isrow(q)
    q=q';
end
if isrow(qd)
    qd=qd';
end

%% Mass, Corriolis, Gravity
g_here    = Tr.g_ini;
J_here    = zeros(6,ndof);
Jd_here   = zeros(6,ndof);
eta_here  = J_here*qd;
dof_start = 1; %starting dof of current piece

for i = 1:N
    
    %Joint
    dof_here = Tr.CVTwists{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    qd_here  = qd(dof_start:dof_start+dof_here-1);
    B_here   = Tr.CVTwists{i}(1).B;
    xi_star  = Tr.CVTwists{i}(1).xi_star;
    
    if dof_here == 0 %fixed joint (N)
        g_joint   = eye(4);
        TgB_here  = zeros(6,ndof);
        TgBd_here = zeros(6,ndof);
    else
        if Tr.VLinks(Tr.LinkIndex(i)).jointtype == 'U' %special case for universal joint. Considered as 2 revolute joints
            % first revolute joint
            xi                                          = B_here(:,1)*q_here(1)+xi_star;
            xid                                         = B_here(:,1)*qd_here(1);
            theta_here                                  = norm(xi(1:3));
            
            thetad_here                                 = (xid(1:3)'*xi(1:3))/theta_here;
            g_joint                                     = joint_expmap(xi);
            
            Tg                                          = variable_Texpmap(1,theta_here,xi);
            TgB_here                                    = zeros(6,ndof);
            TgB_here(:,dof_start)                       = Tg*B_here(:,1);
            
            Tgd                                         = variable_dotTexpmap(1,theta_here,thetad_here,xi,xid);
            TgBd_here                                   = zeros(6,ndof);
            TgBd_here(:,dof_start)                      = dinamico_adj(eta_here)*Tg*B_here(:,1)+Tgd*B_here(:,1);
            
            g_here                                      = g_here*g_joint;
            J_here                                      = dinamico_Adjoint(ginv(g_joint))*(J_here+TgB_here);
            Jd_here                                     = dinamico_Adjoint(ginv(g_joint))*(Jd_here+TgBd_here);
            eta_here                                    = J_here*qd;
            
            % second revolute joint
            xi                                          = B_here(:,2)*q_here(2)+xi_star;
            xid                                         = B_here(:,2)*qd_here(2);
            theta_here                                  = norm(xi(1:3));
            
            thetad_here                                 = (xid(1:3)'*xi(1:3))/theta_here;
            g_joint                                     = joint_expmap(xi);
            
            Tg                                          = variable_Texpmap(1,theta_here,xi);
            TgB_here                                    = zeros(6,ndof);
            TgB_here(:,dof_start+1)                     = Tg*B_here(:,2);
            
            Tgd                                         = variable_dotTexpmap(1,theta_here,thetad_here,xi,xid);
            TgBd_here                                   = zeros(6,ndof);
            TgBd_here(:,dof_start+1)                    = dinamico_adj(eta_here)*Tg*B_here(:,2)+Tgd*B_here(:,2);
        else
            xi                                          = B_here*q_here+xi_star;
            xid                                         = B_here*qd_here;
            theta_here                                  = norm(xi(1:3));
            
            thetad_here                                 = (xid(1:3)'*xi(1:3))/theta_here;
            g_joint                                     = joint_expmap(xi);
            
            Tg                                          = variable_Texpmap(1,theta_here,xi);
            TgB_here                                    = zeros(6,ndof);
            TgB_here(:,dof_start:dof_start+dof_here-1)  = Tg*B_here;
            
            Tgd                                         = variable_dotTexpmap(1,theta_here,thetad_here,xi,xid);
            TgBd_here                                   = zeros(6,ndof);
            TgBd_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*Tg*B_here+Tgd*B_here;
        end
    end
    
    %updating g, Jacobian, Jacobian_dot and eta
    g_here   = g_here*g_joint;
    J_here   = dinamico_Adjoint(ginv(g_joint))*(J_here+TgB_here);
    Jd_here  = dinamico_Adjoint(ginv(g_joint))*(Jd_here+TgBd_here);
    eta_here = J_here*qd;
    
    if Tr.VLinks(Tr.LinkIndex(i)).linktype == 'r'
        
        gi         = Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here     = g_here*gi;
        J_here     = dinamico_Adjoint(ginv(gi))*J_here;
        Jd_here    = dinamico_Adjoint(ginv(gi))*Jd_here;
        eta_here   = J_here*qd;
        
        M_here     = Tr.VLinks(Tr.LinkIndex(i)).Ms;
        C          = C+J_here'*(M_here*Jd_here+dinamico_coadj(eta_here)*M_here*J_here);
        
        % bringing all quantities to the end of rigid link
        gf         = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here     = g_here*gf;
        J_here     = dinamico_Adjoint(ginv(gf))*J_here;
        Jd_here    = dinamico_Adjoint(ginv(gf))*Jd_here;
        eta_here   = J_here*qd;
    end
    
    dof_start = dof_start+dof_here;
    
    for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link
        
        dof_here   = Tr.CVTwists{i}(j+1).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        qd_here    = qd(dof_start:dof_start+dof_here-1);
        Bdof       = Tr.CVTwists{i}(j+1).Bdof;
        Bodr       = Tr.CVTwists{i}(j+1).Bodr;
        xi_star    = Tr.CVTwists{i}(j+1).xi_star;
        gi         = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        lpf        = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Ms         = Tr.VLinks(Tr.LinkIndex(i)).Ms{j};
        Xs         = Tr.VLinks(Tr.LinkIndex(i)).Xs{j};
        Ws         = Tr.VLinks(Tr.LinkIndex(i)).Ws{j};
        nGauss     = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};
        
        q_scale_here  = Tr.q_scale(dof_start:dof_start+dof_here-1);
        doftheta_here = Bdof(1:3)'*(Bodr(1:3)+1);
        q_scale_here(1:doftheta_here) = q_scale_here(1:doftheta_here)*lpf;
        B_scale = repmat(q_scale_here',6*nGauss,1);
        B_Z1 = Tr.CVTwists{i}(j+1).B_Z1./B_scale;
        B_Z2 = Tr.CVTwists{i}(j+1).B_Z2./B_scale;
        
        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here     = g_here*gi;
        J_here     = dinamico_Adjoint(ginv(gi))*J_here;
        Jd_here    = dinamico_Adjoint(ginv(gi))*Jd_here;
        eta_here   = J_here*qd;
        
        for ii = 2:nGauss
            
            H    = (Xs(ii)-Xs(ii-1))*lpf;
            
            B_Z1here      = B_Z1(6*(ii-2)+1:6*(ii-1),:);
            B_Z2here      = B_Z2(6*(ii-2)+1:6*(ii-1),:);
            
            xi_starZ1here = xi_star(6*(ii-2)+1:6*(ii-1),2);
            xi_starZ2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
            
            xi_Z1here = B_Z1here*q_here+xi_starZ1here;
            xi_Z2here = B_Z2here*q_here+xi_starZ2here;
            
            Gamma_here  = (H/2)*(xi_Z1here+xi_Z2here)+...
                          ((sqrt(3)*H^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
            k_here      = Gamma_here(1:3);
            theta_here  = norm(k_here);
            gh          = variable_expmap(theta_here,Gamma_here);
            
            BGamma_here = (H/2)*(B_Z1here+B_Z2here)+...
                          ((sqrt(3)*H^2)/12)*(dinamico_adj(xi_Z1here)*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);
            
            TGamma_here  = variable_Texpmap(1,theta_here,Gamma_here);
            TBGamma_here = zeros(6,ndof);
            TBGamma_here(:,dof_start:dof_start+dof_here-1)  = TGamma_here*BGamma_here;
            
            Gammad_here = BGamma_here*qd_here;
            kd_here     = Gammad_here(1:3);
            
            thetad_here   = (kd_here'*k_here)/theta_here;
            TGammad_here  = variable_dotTexpmap(1,theta_here,thetad_here,Gamma_here,Gammad_here);
            TBGammad_here = zeros(6,ndof);
            TBGammad_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*TGamma_here*BGamma_here+TGammad_here*BGamma_here;
            
            %updating g, Jacobian, Jacobian_dot and eta
            g_here   = g_here*gh;
            J_here   = dinamico_Adjoint(ginv(gh))*(J_here+TBGamma_here);
            Jd_here  = dinamico_Adjoint(ginv(gh))*(Jd_here+TBGammad_here);
            eta_here = J_here*qd;
            
            %integrals evaluation
            if ii<nGauss
                W_here   =Ws(ii);
                Ms_here  =Ms(6*(ii-1)+1:6*ii,:);
                C        =C+lpf*W_here*J_here'*(Ms_here*Jd_here+dinamico_coadj(eta_here)*Ms_here*J_here);
            end
            
        end
        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gf         = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
        g_here     = g_here*gf;
        J_here     = dinamico_Adjoint(ginv(gf))*J_here;
        Jd_here    = dinamico_Adjoint(ginv(gf))*Jd_here;
        eta_here   = J_here*qd;
        
        dof_start  = dof_start+dof_here;
    end
    
end
end

