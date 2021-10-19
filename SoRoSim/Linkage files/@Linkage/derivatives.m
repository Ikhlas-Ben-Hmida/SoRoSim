%Function to calculate qdot and qdotdot at a given time used for ode45
%dynamics
%Last modified by Anup Teejo Mathew - 23/05/2021
function dqdt=derivatives(S,t,qqd,uqt)
  
    t
    K_scale = S.q_scale*S.q_scale';
    %precomputed values
    if S.Damped
        D=S.D.*K_scale;
    else
        D=0;
    end
    K=S.K.*K_scale;
    
    ndof   = S.ndof;
    
    q  = qqd(1:ndof);
    qd = qqd(ndof+1:2*ndof);
    
    n_jact = S.n_jact;
    
    if S.Actuated
        WrenchControlled = S.WrenchControlled;
        i_jactq          = S.i_jactq;
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
    G = S.G;

    N     = S.N;
    n_sig = S.n_sig;
    
 %% Mass, Corriolis, Gravity
    g                 = zeros(4*n_sig,4);
    J                 = zeros(6*n_sig,ndof);
    Jd                = zeros(6*n_sig,ndof);
    eta               = zeros(6*n_sig,1);
    g_here            = S.g_ini;
    J_here            = zeros(6,ndof);
    Jd_here           = zeros(6,ndof);
    eta_here          = J_here*qd;
    dof_start         = 1; %starting dof of current piece
    f                 = 1; %piece number
    i_sig             = 1;

    for i=1:N
        
        %Joint
        dof_here   = S.Vtwists(f).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        qd_here    = qd(dof_start:dof_start+dof_here-1);
        B_here     = S.Vtwists(f).B;
        xi_star    = S.Vtwists(f).xi_star;
        
        if dof_here==0 %fixed joint (N)
            g_joint     = eye(4);
            TgB_here    = zeros(6,ndof);
            TgBd_here   = zeros(6,ndof);
        else
            if S.VLinks(i).jointtype=='U' %special case for universal joint. Considered as 2 revolute joints
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
        
        if S.VLinks(i).linktype=='r'
            
            g0f        = S.g0{f};
            g_here     = g_here*g0f;
            J_here     = dinamico_Adjoint(ginv(g0f))*J_here;
            Jd_here    = dinamico_Adjoint(ginv(g0f))*Jd_here;
            eta_here   = J_here*qd;
            
            g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
            J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
            Jd((i_sig-1)*6+1:i_sig*6,:)   = Jd_here;
            eta((i_sig-1)*6+1:i_sig*6)    = eta_here;
            i_sig                         = i_sig+1;
            
            M_here = S.VLinks(i).Ms;
            M      = M+J_here'*M_here*J_here;
            C      = C+J_here'*(M_here*Jd_here+dinamico_coadj(eta_here)*M_here*J_here);
            if S.Gravity
                F  = F+J_here'*M_here*dinamico_Adjoint(ginv(g_here))*G;
            end
            
            % bringing all quantities to the end of rigid link
            g0f(2:3,4) = -g0f(2:3,4);
            g0f(1,4)   = S.VLinks(i).L-g0f(1,4);
            g_here     = g_here*g0f;
            J_here     = dinamico_Adjoint(ginv(g0f))*J_here;
            Jd_here    = dinamico_Adjoint(ginv(g0f))*Jd_here;
            eta_here   = J_here*qd;
        end
        
        dof_start = dof_start+dof_here;
        f=f+1;

        for j=1:S.VLinks(i).npie-1 %will run only if soft link
            
            dof_here   = S.Vtwists(f).dof;
            q_here     = q(dof_start:dof_start+dof_here-1);
            qd_here    = qd(dof_start:dof_start+dof_here-1);
            B_Z1       = S.Vtwists(f).B_Z1;
            B_Z2       = S.Vtwists(f).B_Z2;
            xi_star    = S.Vtwists(f).xi_star;
            g0f        = S.g0{f};
            lpf        = S.VLinks(i).lp{j};
            Ms         = S.VLinks(i).Ms{j};
            Xs         = S.VLinks(i).Xs{j};
            Ws         = S.VLinks(i).Ws{j};
            nGauss     = S.VLinks(i).nGauss{j};
            
            %updating g, Jacobian, Jacobian_dot and eta at X=0
            g_here        = g_here*g0f;
            J_here        = dinamico_Adjoint(ginv(g0f))*J_here;
            Jd_here       = dinamico_Adjoint(ginv(g0f))*Jd_here;
            eta_here      = J_here*qd;
            
            %scaling of quantities
            Lscale         = lpf;
            g_here(1:3,4)  = g_here(1:3,4)/Lscale;
            J_here(4:6,:)  = J_here(4:6,:)/Lscale;
            Jd_here(4:6,:) = Jd_here(4:6,:)/Lscale;
            eta_here(4:6)  = eta_here(4:6)/Lscale;
            lpf            = 1;
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
                
                H             = (Xs(ii)-Xs(ii-1))*lpf;
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
                    Ms_here(1:3,:)   = Ms_here(1:3,:)/Lscale;
                    Ms_here(4:6,:)   = Ms_here(4:6,:)*Lscale;
                    M        = M+lpf*W_here*J_here'*Ms_here*J_here*Lscale^2; %rescale
                    C        = C+lpf*W_here*J_here'*(Ms_here*Jd_here+dinamico_coadj(eta_here)*Ms_here*J_here)*Lscale^2; %rescale
                    if S.Gravity
                        F    = F+lpf*W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %rescale
                    end
                end
                
            end

            %scaling back quantities
            g_here(1:3,4)  = g_here(1:3,4)*Lscale;
            J_here(4:6,:)  = J_here(4:6,:)*Lscale;
            Jd_here(4:6,:) = Jd_here(4:6,:)*Lscale;
            G              = G*Lscale;

            %updating g, Jacobian, Jacobian_dot and eta at X=L
            g0f(2:3,4)    = -g0f(2:3,4);
            g_here        = g_here*g0f;
            J_here        = dinamico_Adjoint(ginv(g0f))*J_here;
            Jd_here       = dinamico_Adjoint(ginv(g0f))*Jd_here;
            eta_here      = J_here*qd;
            
            dof_start = dof_start+dof_here;
            f=f+1;
        end
        
    end

 %% Point Force
    if S.PointForce
        F = F+ComputePointForce(S,J,t);
    end
    
%% Actuation   
    if S.Actuated
        
        nact = S.nact;
        Bq   = zeros(ndof,nact);
        u    = zeros(nact,1);
        
        %revolute, prismatic, helical joints 
        n_jact           = S.n_jact;
        Bqj1             = S.Bqj1;
        n_1dof           = size(Bqj1,2);
        Bq(:,1:n_1dof)   = S.Bqj1;
        
        %for other joints
        i_jact           = S.i_jact;
        i_u              = n_1dof+1;
        i_jactq          = S.i_jactq;
        
        
        n_ljact=length(i_jact);
        for iii=i_u:n_ljact
            i=i_jact(iii);
            if S.VLinks(i).jointtype=='U'
                Bq(i_jactq(i_u:i_u+1),i_u:i_u+1) = [1 0;0 1];
                i_u = i_u+2;
            elseif S.VLinks(i).jointtype=='C'
                Bq(i_jactq(i_u:i_u+1),i_u:i_u+1) = [1 0;0 1];
                i_u = i_u+2;
            elseif S.VLinks(i).jointtype=='A'
                i_sig = 1;
                f     = 1;
                for ii=1:i-1
                    i_sig = i_sig+1;
                    for jj=1:S.VLinks(ii).npie-1
                        i_sig = i_sig+1+S.VLinks(ii).nGauss{jj};
                    end
                    if S.VLinks(ii).linktype=='r'
                        i_sig = i_sig+1;
                    end
                    f=f+S.VLinks(ii).npie;
                end
                J_here = J((i_sig-1)*6+1:i_sig*6,:);
                S_here = J_here(:,i_jactq(i_u:i_u+2));
                B_here = S.Vtwists(f).B;
                
                Bq(i_jactq(i_u:i_u+2),i_u:i_u+2) = S_here'*B_here;
                i_u = i_u+3;
            elseif S.VLinks(i).jointtype=='S'
                i_sig = 1;
                for ii=1:i-1
                    i_sig = i_sig+1;
                    for jj=1:S.VLinks(ii).npie-1
                        i_sig = i_sig+1+S.VLinks(ii).nGauss{jj};
                    end
                    if S.VLinks(ii).linktype=='r'
                        i_sig = i_sig+1;
                    end
                end
                J_here = J((i_sig-1)*6+1:i_sig*6,:);
                S_here = J_here(:,i_jactq(i_u:i_u+2));
                B_here = [eye(3);zeros(3,3)];
                
                Bq(i_jactq(i_u:i_u+2),i_u:i_u+2) = S_here'*B_here;
                i_u = i_u+3;
            else %free joint
                i_sig = 1;
                for ii=1:i-1
                    i_sig = i_sig+1;
                    for jj=1:S.VLinks(ii).npie-1
                        i_sig = i_sig+1+S.VLinks(ii).nGauss{jj};
                    end
                    if S.VLinks(ii).linktype=='r'
                        i_sig = i_sig+1;
                    end
                end
                J_here = J((i_sig-1)*6+1:i_sig*6,:);
                S_here = J_here(:,i_jactq(i_u:i_u+5));
                B_here = eye(6);
                
                Bq(i_jactq(i_u:i_u+5),i_u:i_u+5) = S_here'*B_here;
                i_u = i_u+6;
            end    
        end

        %Cable actuation
        n_sact=S.n_sact;
        
        for ii=1:n_sact
            
            dcii = cell(1,N); dcpii = cell(1,N); Sdivii = cell(1,N); Edivii = cell(1,N);
            
            for i=1:N
                dcii{i}   = S.dc{ii,i};
                dcpii{i}  = S.dcp{ii,i};
                Sdivii{i} = S.Sdiv{ii,i};
                Edivii{i} = S.Ediv{ii,i};
            end
            Insideii     = S.Inside{ii};
            
            if Insideii
                Bq(:,n_jact+ii) = ComputeCableActuation(S,dcii,dcpii,Sdivii,Edivii,q);
            else
                Bq(:,n_jact+ii) = ComputeCableActuation2(S,dcii,Sdivii,Edivii,J,g);
            end

        end
        
        if ~S.CAS
            WrenchControlled = S.WrenchControlled;
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

%% Custom External Force, Actuation or Actuator Strength  
    if S.CEFP||S.CAP||S.CAS

        J_scale  = repmat(S.q_scale',6,1);
        i_sig    = 1;
        g_act    = zeros(n_sig*4,4);
        J_act    = zeros(n_sig*6,ndof);
        eta_act  = zeros(n_sig*6,1);
        Jd_act   = zeros(n_sig*6,ndof);
        for i=1:N
            % joint
            g_act((i_sig-1)*4+1:i_sig*4,:)  = g((i_sig-1)*4+1:i_sig*4,:);
            J_act((i_sig-1)*6+1:i_sig*6,:)  = J((i_sig-1)*6+1:i_sig*6,:).*J_scale;
            eta_act((i_sig-1)*6+1:i_sig*6)  = eta((i_sig-1)*6+1:i_sig*6);
            Jd_act((i_sig-1)*6+1:i_sig*6,:) = Jd((i_sig-1)*6+1:i_sig*6,:).*J_scale;
            i_sig = i_sig+1;
            if S.VLinks(i).linktype=='r'
                g_act((i_sig-1)*4+1:i_sig*4,:)  = g((i_sig-1)*4+1:i_sig*4,:);
                J_act((i_sig-1)*6+1:i_sig*6,:)  = J((i_sig-1)*6+1:i_sig*6,:).*J_scale;
                eta_act((i_sig-1)*6+1:i_sig*6)  = eta((i_sig-1)*6+1:i_sig*6);
                Jd_act((i_sig-1)*6+1:i_sig*6,:) = Jd((i_sig-1)*6+1:i_sig*6,:).*J_scale;
                i_sig = i_sig+1;
            end
            for j=1:S.VLinks(i).npie-1
                Lscale = S.VLinks(i).lp{j};
                for ii=1:S.VLinks(i).nGauss{j}

                g_here        = g((i_sig-1)*4+1:i_sig*4,:);
                g_here(1:3,4) = g_here(1:3,4)*Lscale;
                g_act((i_sig-1)*4+1:i_sig*4,:)  = g_here;

                J_here        = J((i_sig-1)*6+1:i_sig*6,:);
                J_here(4:6,:) = J_here(4:6,:)*Lscale;
                J_act((i_sig-1)*6+1:i_sig*6,:)  = J_here.*J_scale;

                eta_here      = eta((i_sig-1)*6+1:i_sig*6);
                eta_here(4:6) = eta_here(4:6)*Lscale;
                eta_act((i_sig-1)*6+1:i_sig*6)  = eta_here;

                Jd_here        = Jd((i_sig-1)*6+1:i_sig*6,:);
                Jd_here(4:6,:) = Jd_here(4:6,:)*Lscale;
                Jd_act((i_sig-1)*6+1:i_sig*6,:)  = Jd_here.*J_scale;
                i_sig = i_sig+1;
                end
            end
        end
        if S.CEFP
            Fext = CustomExtForce(S,q.*S.q_scale,g_act,J_act,t,qd.*S.q_scale,eta_act,Jd_act);
            [F_custom,M_custom,C_custom] = CustomExtForce_Qspace(S,J,Fext,Jd,eta);
            M = M+M_custom;
            F = F+F_custom;
            C = C+C_custom;
        end
        if S.CAP
            Fact = CustomActuation(S,q.*S.q_scale,g_act,J_act,t,qd.*S.q_scale,eta_act,Jd_act);
            TauC = CustomActuation_Qspace(S,Fact);
            F    = F+TauC; %added with F
        end
        if S.CAS
            u = CustomActuatorStrength(S,q.*S.q_scale,g_act,J_act,t,qd.*S.q_scale,eta_act,Jd_act);
        end
        
    end
%%    

    qdd = M\(Bq*u+F-K*q-(C+D)*qd);
    if S.Actuated
        for i=1:n_jact
            if ~WrenchControlled(i)
                qdd(i_jactq(i)) = u(i);
            end
        end
    end
    
    dqdt=[qd;qdd];
end