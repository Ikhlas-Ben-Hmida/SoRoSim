%Computes the generalized external force acting on the linkage in Q space
%Last modified by Anup Teejo Mathew - 24/05/2021
function F = GeneralizedExternalForce(Tr,q)

    if isrow(q)
        q=q';
    end

    N      = Tr.N;
    ndof   = Tr.ndof;
    nsig   = Tr.nsig;
    
    
    i_sig        = 1;
    J            = zeros(6*nsig,ndof);
    g            = zeros(4*nsig,4);
    g_here       = Tr.g_ini;
    J_here       = zeros(6,ndof);
    dof_start    = 1; %starting dof of current piece
    F            = zeros(ndof,1);%external force 
    G            = Tr.G;
    
    for i=1:N

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
        g_here                        = g_here*g_joint;
        J_here                        = dinamico_Adjoint(ginv(g_joint))*(J_here+TgB_here);
        g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
        J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
        i_sig                         = i_sig+1;
        if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
            
            gi         = Tr.VLinks(Tr.LinkIndex(i)).gi;
            g_here     = g_here*gi;
            J_here     = dinamico_Adjoint(ginv(gi))*J_here;
            
            g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
            J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
            i_sig                         = i_sig+1;
                        
            M_here=Tr.VLinks(Tr.LinkIndex(i)).Ms;
            if Tr.Gravity
                F = F+J_here'*M_here*dinamico_Adjoint(ginv(g_here))*G;
            end
            % bringing all quantities to the end of rigid link
            gf         = Tr.VLinks(Tr.LinkIndex(i)).gf;
            g_here     = g_here*gf;
            J_here     = dinamico_Adjoint(ginv(gf))*J_here;
        end

        dof_start=dof_start+dof_here;

        for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link

            dof_here   = Tr.CVTwists{i}(j+1).dof;
            q_here     = q(dof_start:dof_start+dof_here-1);
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
            g_here        = g_here*gi;
            J_here        = dinamico_Adjoint(ginv(gi))*J_here;
            
            g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
            J((i_sig-1)*6+1:i_sig*6,:)    = J_here;
            i_sig                         = i_sig+1;
            
            for ii=2:nGauss

                H    = (Xs(ii)-Xs(ii-1))*lpf;

                B_Z1here      = B_Z1(6*(ii-2)+1:6*(ii-1),:);
                B_Z2here      = B_Z2(6*(ii-2)+1:6*(ii-1),:);
                
                xi_starZ1here = xi_star(6*(ii-2)+1:6*(ii-1),2);
                xi_starZ2here = xi_star(6*(ii-2)+1:6*(ii-1),3);

                xi_Z1here     = B_Z1here*q_here+xi_starZ1here;
                xi_Z2here     = B_Z2here*q_here+xi_starZ2here;

                Gamma_here    = (H/2)*(xi_Z1here+xi_Z2here)+...
                                ((sqrt(3)*H^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
                k_here        = Gamma_here(1:3);
                theta_here    = norm(k_here);
                gh            = variable_expmap(theta_here,Gamma_here);

                BGamma_here       = (H/2)*(B_Z1here+B_Z2here)+...
                                    ((sqrt(3)*H^2)/12)*(dinamico_adj(xi_Z1here)*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);

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
                    if Tr.Gravity
                        F = F+lpf*W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G;
                    end
                end

            end
            %updating g, Jacobian, Jacobian_dot and eta at X=L
            gf            = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
            g_here        = g_here*gf;
            J_here        = dinamico_Adjoint(ginv(gf))*J_here;
%             J_here        = ginv(gf)*J_here;


            dof_start = dof_start+dof_here;
        end

    end
    %% Point Force
    if Tr.PointForce
        F = F+ComputePointForce(Tr,J,g);
    end

end

