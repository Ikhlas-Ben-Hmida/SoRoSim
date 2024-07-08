%Function that calculates the forward kinematics of the linkage from the base to the linkage tip at every significant points(24.05.2021)
%in SI units Last modified by Anup Teejo Mathew 02.03.2022
function g = FwdKinematics(Tr,q,i_here,j_here) %i_here is link index, j_here is division (0 for joint)

full = false;
if nargin==4
    if j_here==0
        nsig = 1;
    else
        nsig = Tr.CVTwists{i_here}(j_here+1).nip;
    end
elseif nargin==3
    if Tr.VLinks(Tr.LinkIndex(i_here)).linktype=='s'
        nsig = 1;
        for j=1:Tr.VLinks(Tr.LinkIndex(i_here)).npie-1
            nsig = nsig+Tr.CVTwists{i_here}(j+1).nip;
        end
    else
        nsig = 2; %joint and CM
    end
else
    full  = true;
    nsig  = Tr.nsig;
end
g     = zeros(4*nsig,4); 
i_sig = 1;

if isrow(q)
    q=q';
end

N         = Tr.N;
g_ini     = Tr.g_ini;
g_Ltip    = repmat(eye(4),N,1);
iLpre     = Tr.iLpre;

dof_start = 1;                         %starting dof of current piece


for i = 1:N
    
    if Tr.iLpre(i)>0
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
        xi      = B_here*q_here+xi_star;
        g_joint = variable_expmap_g(xi);
    end

    g_here                     = g_here*g_joint;
    if full||(i==i_here&&nargin==3)||(i==i_here&&j_here==0)
        g((i_sig-1)*4+1:i_sig*4,:) = g_here;
        i_sig                      = i_sig+1;
    end
    
    if Tr.VLinks(Tr.LinkIndex(i)).linktype == 'r'
        
        gi                         = Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here                     = g_here*gi;
        if full||(i==i_here&&nargin==3)
            g((i_sig-1)*4+1:i_sig*4,:) = g_here;
            i_sig                      = i_sig+1;
        end
        % bringing all quantities to the end of rigid link
        gf     = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here = g_here*gf;
    end
    
    dof_start = dof_start+dof_here;
    
    for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        
        dof_here = Tr.CVTwists{i}(j+1).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);
        xi_star  = Tr.CVTwists{i}(j+1).xi_star;
        gi       = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        ld      = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Xs       = Tr.CVTwists{i}(j+1).Xs;
        nip      = Tr.CVTwists{i}(j+1).nip;
        
        if Tr.Z_order==4
            B_Z1       = Tr.CVTwists{i}(j+1).B_Z1;
            B_Z2       = Tr.CVTwists{i}(j+1).B_Z2;
        else %order 2
            B_Z        = Tr.CVTwists{i}(j+1).B_Z;
        end
        
        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here = g_here*gi;
        if full||(i==i_here&&nargin==3)||(i==i_here&&j==j_here)
            g((i_sig-1)*4+1:i_sig*4,:) = g_here;
            i_sig                      = i_sig+1;
        end
        Lscale = ld; 
        
        for ii = 2:nip

            H    = Xs(ii)-Xs(ii-1);
            
            if Tr.Z_order==4
                
                xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1),2);
                xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
                xi_Z1here(1:3) = xi_Z1here(1:3)*Lscale; %scaling
                xi_Z2here(1:3) = xi_Z2here(1:3)*Lscale;
                
                if dof_here>0
                    if dof_here>0
                        B_Z1here  = B_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
                        B_Z2here  = B_Z2(6*(ii-2)+1:6*(ii-1),:);
                        xi_Z1here = B_Z1here*q_here+xi_Z1here;
                        xi_Z2here = B_Z2here*q_here+xi_Z2here;
                    end
                end
                ad_xi_Z1here = dinamico_adj(xi_Z1here);  
                Gamma_here   = (H/2)*(xi_Z1here+xi_Z2here)+...
                                ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;
                      
            else % order 2
                xi_Zhere = xi_star(6*(ii-2)+1:6*(ii-1),4);
                xi_Zhere(1:3) = xi_Zhere(1:3)*Lscale; %scaling
                
                if dof_here>0
                    B_Zhere  = B_Z(6*(ii-2)+1:6*(ii-1),:);%note this step
                    if dof_here>0
                        xi_Zhere = B_Zhere*q_here+xi_Zhere;
                    end
                end
                
                Gamma_here  = H*xi_Zhere;

            end
            Gamma_here(4:6) = Gamma_here(4:6)*Lscale;
            gh              = variable_expmap_g(Gamma_here);

            %updating g, Jacobian, Jacobian_dot and eta
            
            g_here                     = g_here*gh;
            if full||(i==i_here&&nargin==3)||(i==i_here&&j==j_here)
                g((i_sig-1)*4+1:i_sig*4,:) = g_here;
                i_sig                      = i_sig+1;
            end

        end

        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
        g_here    = g_here*gf;
        dof_start = dof_start+dof_here;
        
    end
    g_Ltip((i-1)*4+1:i*4,:) = g_here; 
end

end


