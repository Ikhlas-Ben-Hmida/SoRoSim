%Function that calculates the forward kinematics of the linkage from the base to the linkage tip at every significant points(24.05.2021)
%Last modified by Anup Teejo Mathew 02.03.2022
function g = FwdKinematics(Tr,q)

if isrow(q)
    q=q';
end

N         = Tr.N;
nsig      = Tr.nsig;
g_ini     = Tr.g_ini;
g_Ltip    = repmat(eye(4),N,1);
iLpre     = Tr.iLpre;
LinkIndex = Tr.LinkIndex;

g         = zeros(4*nsig,4);
dof_start = 1;                         %starting dof of current piece
i_sig     = 1;

for i = 1:N
    
    if iLpre(i)>0
        g_here=g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
    else
        g_here=g_ini((i-1)*4+1:i*4,:);
    end
    
    %Joint
    dof_here = Tr.CVTwists{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    B_here   = Tr.CVTwists{i}(1).B;
    xi_star  = Tr.CVTwists{i}(1).xi_star;
    
    if dof_here == 0                   %fixed joint (N)
        g_joint  = eye(4);
    else
        if Tr.VLinks(LinkIndex(i)).jointtype=='U'  %special case for universal joint. Considered as 2 revolute joints
            % first revolute joint
            xi                                         = B_here(:,1)*q_here(1)+xi_star;
            g_joint                                    = joint_expmap(xi);
            g_here                                     = g_here*g_joint;
            % second revolute joint
            xi                                         = B_here(:,2)*q_here(2)+xi_star;
            g_joint                                    = joint_expmap(xi);
        else
            xi                                         = B_here*q_here+xi_star;
            g_joint                                    = joint_expmap(xi);
        end
    end

    g_here                     = g_here*g_joint;
    g((i_sig-1)*4+1:i_sig*4,:) = g_here;
    i_sig                      = i_sig+1;
    
    if Tr.VLinks(LinkIndex(i)).linktype == 'r'
        
        gi                         = Tr.VLinks(LinkIndex(i)).gi;
        g_here                     = g_here*gi;
        g((i_sig-1)*4+1:i_sig*4,:) = g_here;
        i_sig                      = i_sig+1;
        % bringing all quantities to the end of rigid link
        gf     = Tr.VLinks(LinkIndex(i)).gf;
        g_here = g_here*gf;
    end
    
    dof_start = dof_start+dof_here;
    
    for j = 1:Tr.VLinks(LinkIndex(i)).npie-1
        
        dof_here = Tr.CVTwists{i}(j+1).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);
        xi_star  = Tr.CVTwists{i}(j+1).xi_star;
        gi       = Tr.VLinks(LinkIndex(i)).gi{j};
        Bdof     = Tr.CVTwists{i}(j+1).Bdof;
        Bodr     = Tr.CVTwists{i}(j+1).Bodr;
        lpf      = Tr.VLinks(LinkIndex(i)).lp{j};
        Xs       = Tr.VLinks(LinkIndex(i)).Xs{j};
        Z1       = 1/2-sqrt(3)/6;      %Zanna quadrature coefficient
        Z2       = 1/2+sqrt(3)/6;      %Zanna quadrature coefficient
        B_Z1here = zeros(6,dof_here);
        B_Z2here = zeros(6,dof_here);
        nGauss   = Tr.VLinks(LinkIndex(i)).nGauss{j};
        
        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here                     = g_here*gi;
        g((i_sig-1)*4+1:i_sig*4,:) = g_here;
        i_sig                      = i_sig+1;
        
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
            
            xi_Z1here                  = B_Z1here*q_here+xi_starZ1here;
            xi_Z2here                  = B_Z2here*q_here+xi_starZ2here;
 
            Gamma_here                 = (H/2)*(xi_Z1here+xi_Z2here)+...
                                         ((sqrt(3)*H^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
            k_here                     = Gamma_here(1:3);
            theta_here                 = norm(k_here);
            gh                         = variable_expmap(theta_here,Gamma_here);

            %updating g, Jacobian, Jacobian_dot and eta
            g_here                     = g_here*gh;
            g((i_sig-1)*4+1:i_sig*4,:) = g_here;
            i_sig                      = i_sig+1;
            
        end

        gf        = Tr.VLinks(LinkIndex(i)).gf{j};
        g_here    = g_here*gf;
        dof_start = dof_start+dof_here;
    end
    g_Ltip((i-1)*4+1:i*4,:) = g_here; 
end

end


