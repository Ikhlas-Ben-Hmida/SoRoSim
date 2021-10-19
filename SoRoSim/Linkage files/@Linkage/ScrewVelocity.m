%Function that calculates the screw velocity (eta) at every significant point
%(24.05.2021)

function eta = ScrewVelocity(S,q,qd)


if isrow(q)
    q=q';
end
if isrow(qd)
    qd=qd';
end



N         = S.N;
n_sig     = S.n_sig;
ndof      = S.ndof;

%% Mass, Corriolis, Gravity
eta       = zeros(6*n_sig,1);
g_here    = S.g_ini;
J_here    = zeros(6,ndof);
Jd_here   = zeros(6,ndof);
eta_here  = J_here*qd;
dof_start = 1; %starting dof of current piece
f         = 1; %piece number
i_sig     = 1;

for i = 1:N
    %Joint
    dof_here = S.Vtwists(f).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    qd_here  = qd(dof_start:dof_start+dof_here-1);
    B_here   = S.Vtwists(f).B;
    xi_star  = S.Vtwists(f).xi_star;
    
    if dof_here == 0 %fixed joint (N)
        g_joint   = eye(4);
        TgB_here  = zeros(6,ndof);
        TgBd_here = zeros(6,ndof);
    else
        if S.VLinks(i).jointtype == 'U' %special case for universal joint. Considered as 2 revolute joints
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
    g_here                       = g_here*g_joint;
    J_here                       = dinamico_Adjoint(ginv(g_joint))*(J_here+TgB_here);
    Jd_here                      = dinamico_Adjoint(ginv(g_joint))*(Jd_here+TgBd_here);
    eta_here                     = J_here*qd;
    
    eta((i_sig-1)*6+1:i_sig*6) = eta_here;
    i_sig                        = i_sig+1;
    
    if S.VLinks(i).linktype == 'r'
        g0f                          = S.g0{f};
        g_here                       = g_here*g0f;
        J_here                       = dinamico_Adjoint(ginv(g0f))*J_here;
        Jd_here                      = dinamico_Adjoint(ginv(g0f))*Jd_here;
        eta_here                     = J_here*qd;
        
        eta((i_sig-1)*6+1:i_sig*6,:) = eta_here;
        i_sig                        = i_sig+1;
        
        % bringing all quantities to the end of rigid link
        g0f(2:3,4)                   = -g0f(2:3,4);
        g0f(1,4)                     = S.VLinks(i).L-g0f(1,4);
        g_here                       = g_here*g0f;
        J_here                       = dinamico_Adjoint(ginv(g0f))*J_here;
        Jd_here                      = dinamico_Adjoint(ginv(g0f))*Jd_here;
        eta_here                     = J_here*qd;
    end
    
    dof_start = dof_start+dof_here;
    f         = f+1;
    
    for j = 1:S.VLinks(i).npie-1 %will run only if soft link
        dof_here = S.Vtwists(f).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);
        qd_here  = qd(dof_start:dof_start+dof_here-1);
        Bdof     = S.Vtwists(f).Bdof;
        Bodr     = S.Vtwists(f).Bodr;
        xi_star  = S.Vtwists(f).xi_star;
        g0f      = S.g0{f};
        lpf      = S.VLinks(i).lp{j};
        Xs       = S.VLinks(i).Xs{j};
        Z1       = 1/2-sqrt(3)/6;      %Zanna quadrature coefficient
        Z2       = 1/2+sqrt(3)/6;      %Zanna quadrature coefficient
        B_Z1here = zeros(6,dof_here);
        B_Z2here = zeros(6,dof_here);
        nGauss                       = S.VLinks(i).nGauss{j};
        
        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here                       = g_here*g0f;
        J_here                       = dinamico_Adjoint(ginv(g0f))*J_here;
        Jd_here                      = dinamico_Adjoint(ginv(g0f))*Jd_here;
        eta_here                     = J_here*qd;
        
        eta((i_sig-1)*6+1:i_sig*6,:) = eta_here;
        i_sig                        = i_sig+1;
        
        for ii = 2:nGauss

            x    = Xs(ii-1)*lpf;
            H    = (Xs(ii)-Xs(ii-1))*lpf;
            x_Z1 = x+Z1*H;
            x_Z2 = x+Z2*H;

            for jj = 1:6
                for k = 1:Bdof(jj)*Bodr(jj)+Bdof(jj)
                    kk              = sum(Bdof(1:jj-1).*Bodr(1:jj-1))+sum(Bdof(1:jj-1))+k;
                    B_Z1here(jj,kk) = x_Z1^(k-1);
                    B_Z2here(jj,kk) = x_Z2^(k-1);
                end
            end
            
            xi_starZ1here = xi_star(6*(ii-2)+1:6*(ii-1),2);
            xi_starZ2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
            
            xi_Z1here                                       = B_Z1here*q_here+xi_starZ1here;
            xi_Z2here                                       = B_Z2here*q_here+xi_starZ2here;
            
            Gamma_here                                      = (H/2)*(xi_Z1here+xi_Z2here)+...
                ((sqrt(3)*H^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
            
            k_here                                          = Gamma_here(1:3);
            theta_here                                      = norm(k_here);
            gh                                              = variable_expmap(theta_here,Gamma_here);
            
            BGamma_here                                     = (H/2)*(B_Z1here+B_Z2here)+...
                ((sqrt(3)*H^2)/12)*(dinamico_adj(xi_Z1here)*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);
            
            TGamma_here                                     = variable_Texpmap(1,theta_here,Gamma_here);
            TBGamma_here                                    = zeros(6,ndof);
            TBGamma_here(:,dof_start:dof_start+dof_here-1)  = TGamma_here*BGamma_here;
            
            Gammad_here                                     = BGamma_here*qd_here;
            kd_here                                         = Gammad_here(1:3);
            
            thetad_here                                     = (kd_here'*k_here)/theta_here;
            TGammad_here                                    = variable_dotTexpmap(1,theta_here,thetad_here,Gamma_here,Gammad_here);
            TBGammad_here                                   = zeros(6,ndof);
            TBGammad_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*TGamma_here*BGamma_here+TGammad_here*BGamma_here;
            
            %updating g, Jacobian, Jacobian_dot and eta
            g_here                       = g_here*gh;
            J_here                       = dinamico_Adjoint(ginv(gh))*(J_here+TBGamma_here);
            Jd_here                      = dinamico_Adjoint(ginv(gh))*(Jd_here+TBGammad_here);
            eta_here                     = J_here*qd;
            
            eta((i_sig-1)*6+1:i_sig*6,:) = eta_here;
            i_sig                        = i_sig+1;
            
        end
        %updating g, Jacobian, Jacobian_dot and eta at X=L
        g0f(2:3,4) = -g0f(2:3,4);
        g_here     = g_here*g0f;
        J_here     = dinamico_Adjoint(ginv(g0f))*J_here;
        Jd_here    = dinamico_Adjoint(ginv(g0f))*Jd_here;
        eta_here   = J_here*qd;
        
        dof_start  = dof_start+dof_here;
        f          = f+1;
    end
    
end
end

