%Evaluates the static equilibrium equation given by: K(q-q*)-Tau-F=0 for fsolve
%Last modified by Anup Teejo Mathew 23/05/2021
function E=Equilibrium(S,qu,uq)

    N      = S.N;
    ndof   = S.ndof;
    n_sig  = S.n_sig;
    G      = S.G;
    q      = qu;
    
    if S.Actuated
        n_jact           = S.n_jact;
        WrenchControlled = S.WrenchControlled;
        i_jactq          = S.i_jactq;
        for i=1:n_jact
            if ~WrenchControlled(i)
                q(i_jactq(i))   = uq(i);
            end
        end
    end
    K_scale    = S.q_scale*S.q_scale';
    K          = S.K.*K_scale;
    J          = zeros(6*n_sig,ndof);
    g          = zeros(4*n_sig,4);
    f          = 1; 
    g_here     = S.g_ini;
    J_here     = zeros(6,ndof);
    dof_start  = 1; %starting dof of current piece
    i_sig      = 1;
    F          = zeros(ndof,1);%external force
    
    for i=1:N

        %Joint
        dof_here   = S.Vtwists(f).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        B_here     = S.Vtwists(f).B;
        xi_star    = S.Vtwists(f).xi_star;

        if dof_here==0 %fixed joint (N)
            g_joint     = eye(4);
            TgB_here    = zeros(6,ndof);
        else
            if S.VLinks(i).jointtype=='U' %special case for universal joint. Considered as 2 revolute joints
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
        
        if S.VLinks(i).linktype=='r'
            
            g0f        = S.g0{f};
            g_here     = g_here*g0f;
            J_here     = dinamico_Adjoint(ginv(g0f))*J_here;
            
            g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
            J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
            i_sig                         = i_sig+1;
            
            M_here = S.VLinks(i).Ms;
            if S.Gravity
                F = F+J_here'*M_here*dinamico_Adjoint(ginv(g_here))*G;
            end
            % bringing all quantities to the end of rigid link
            g0f(2:3,4) = -g0f(2:3,4);
            g0f(1,4)   = S.VLinks(i).L-g0f(1,4);
            g_here     = g_here*g0f;
            J_here     = dinamico_Adjoint(ginv(g0f))*J_here;
        end

        dof_start = dof_start+dof_here;
        f=f+1;

        for j=1:S.VLinks(i).npie-1 %will run only if soft link

            dof_here   = S.Vtwists(f).dof;
            q_here     = q(dof_start:dof_start+dof_here-1);
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

            %scaling of quantities
            Lscale        = lpf;
            g_here(1:3,4) = g_here(1:3,4)/Lscale;
            J_here(4:6,:) = J_here(4:6,:)/Lscale;
            lpf           = 1;
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

                h             = (Xs(ii)-Xs(ii-1))*lpf;
                Gamma_here    = (h/2)*(xi_Z1here+xi_Z2here)+...
                                ((sqrt(3)*h^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
                k_here        = Gamma_here(1:3);
                theta_here    = norm(k_here);
                gh            = variable_expmap(theta_here,Gamma_here);

                BGamma_here       = (h/2)*(B_Z1here+B_Z2here)+...
                                   ((sqrt(3)*h^2)/12)*(dinamico_adj(xi_Z1here)*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);

                TGamma_here                                    = variable_Texpmap(1,theta_here,Gamma_here);
                TBGamma_here                                   = zeros(6,ndof);
                TBGamma_here(:,dof_start:dof_start+dof_here-1) = TGamma_here*BGamma_here;

                %updating g, Jacobian, Jacobian_dot and eta
                g_here     = g_here*gh;
                J_here     = dinamico_Adjoint(ginv(gh))*(J_here+TBGamma_here);

                g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
                J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
                i_sig                         = i_sig+1;
                %integrals evaluation
                if ii<nGauss
                    W_here  = Ws(ii);
                    Ms_here = Ms(6*(ii-1)+1:6*ii,:);
                    %scaling
                    Ms_here(1:3,:)   = Ms_here(1:3,:)/Lscale;
                    Ms_here(4:6,:)   = Ms_here(4:6,:)*Lscale;
                    if S.Gravity
                        F = F+(lpf)*W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %scaled back for addition
                    end
                end

            end

            %scaling back quantities
            g_here(1:3,4) = g_here(1:3,4)*Lscale;
            J_here(4:6,:) = J_here(4:6,:)*Lscale;
            G             = G*Lscale;

            %updating g, Jacobian, Jacobian_dot and eta at X=L
            g0f(2:3,4)    = -g0f(2:3,4);
            g_here        = g_here*g0f;
            J_here        = dinamico_Adjoint(ginv(g0f))*J_here;

            dof_start = dof_start+dof_here;
            f=f+1;
        end

    end
    %% Point Force
    if S.PointForce
        F = F+ComputePointForce(S,J,0);
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
        i_tau            = n_1dof+1;
        i_jactq          = S.i_jactq;
        
        for i=i_jact(i_tau:end)
            if S.VLinks(i).jointtype=='U'
                Bq(i_jactq(i_tau:i_tau+1),i_tau:i_tau+1) = [1 0;0 1];
                i_tau = i_tau+2;
            elseif S.VLinks(i).jointtype=='C'
                Bq(i_jactq(i_tau:i_tau+1),i_tau:i_tau+1) = [1 0;0 1];
                i_tau = i_tau+2;
            elseif S.VLinks(i).jointtype=='A'
                i_sig=1;
                f    =1;
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
                S_here = J_here(:,i_jactq(i_tau:i_tau+2));
                B_here = S.Vtwists(f).B;
                
                Bq(i_jactq(i_tau:i_tau+2),i_tau:i_tau+2) = S_here'*B_here;
                i_tau = i_tau+3;
            elseif S.VLinks(i).jointtype=='S'
                i_sig=1;
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
                S_here = J_here(:,i_jactq(i_tau:i_tau+2));
                B_here = [eye(3);zeros(3,3)];
                
                Bq(i_jactq(i_tau:i_tau+2),i_tau:i_tau+2) = S_here'*B_here;
                i_tau = i_tau+3;
            else %free joint
                i_sig=1;
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
                S_here = J_here(:,i_jactq(i_tau:i_tau+5));
                B_here = eye(6);
                
                Bq(i_jactq(i_tau:i_tau+5),i_tau:i_tau+5) = S_here'*B_here;
                i_tau = i_tau+6;
            end    
        end
             
        %cable actuation
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
                    u(i)           =qu(i_jactq(i));
                else
                    u(i)           =uq(i);
                end
            end

            for i=n_jact+1:nact
                u(i)=uq(i);
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

    for i=1:N
        % joint
        g_act((i_sig-1)*4+1:i_sig*4,:)  = g((i_sig-1)*4+1:i_sig*4,:);
        J_act((i_sig-1)*6+1:i_sig*6,:)  = J((i_sig-1)*6+1:i_sig*6,:).*J_scale;
        i_sig = i_sig+1;
        if S.VLinks(i).linktype=='r'
            g_act((i_sig-1)*4+1:i_sig*4,:)  = g((i_sig-1)*4+1:i_sig*4,:);
            J_act((i_sig-1)*6+1:i_sig*6,:)  = J((i_sig-1)*6+1:i_sig*6,:).*J_scale;
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

            i_sig = i_sig+1;
            end
        end
    end
    
    if S.CEFP
        Fext = CustomExtForce(S,q.*S.q_scale,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof));
        [F_custom,~,~] = CustomExtForce_Qspace(S,J,Fext,zeros(6*n_sig,ndof),zeros(6*n_sig,1));
        F              = F+F_custom;
    end
    if S.CAP
        Fact = CustomActuation(S,q.*S.q_scale,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof));
        TauC = CustomActuation_Qspace(S,Fact);
        F    = F+TauC; %added with F
    end
    if S.CAS
        u = CustomActuatorStrength(S,q.*S.q_scale,g_act,J_act,0,zeros(ndof,1),zeros(6*n_sig,1),zeros(6*n_sig,ndof));
    end
end
    %% Equilibrium
    E=K*q-Bq*u-F;
    E=E*1e5;

end