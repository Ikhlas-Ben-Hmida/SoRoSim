%Function that calculates the jacobian at every significant points
% in SI units (Last modified 18.04.2023, Anup)
function J = Jacobian(Tr,q,i_here,j_here)%i_here is link index, j_here is division (0 for joint)

if isrow(q)
    q=q';
end


full = false;
if nargin==4
    if j_here==0
        nsig = 1;
    else
        nsig = Tr.CVTwists{i_here}(j_here+1).nip; %j_here>1 is allowed only for soft links
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
ndof  = Tr.ndof;
J     = zeros(6*nsig,ndof);
i_sig = 1;

N         = Tr.N;

g_ini     = Tr.g_ini;
g_Ltip    = repmat(eye(4),N,1);
J_Ltip    = repmat(zeros(6,ndof),N,1);
iLpre     = Tr.iLpre;

dof_start = 1; %starting dof of current piece

for i = 1:N
    
    if iLpre(i)>0
        g_here       = g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        J_here       = Ad_g_ini_inv*J_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
    else
        g_here   = g_ini((i-1)*4+1:i*4,:);
        J_here   = zeros(6,ndof);
    end
    
    %Joint
    dof_here = Tr.CVTwists{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    B_here   = Tr.CVTwists{i}(1).B;
    xi_star  = Tr.CVTwists{i}(1).xi_star;
    
    if dof_here == 0 %fixed joint (N)
        g_joint  = eye(4);
        TgB_here = zeros(6,ndof);
    else
        xi           = B_here*q_here+xi_star;
        [g_joint,Tg] = variable_expmap_gTg(xi);

        TgB_here                                   = zeros(6,ndof);
        TgB_here(:,dof_start:dof_start+dof_here-1) = Tg*B_here;
    end
    
    %updating g, Jacobian, Jacobian_dot and eta
    g_here = g_here*g_joint;
    J_here = dinamico_Adjoint(ginv(g_joint))*(J_here+TgB_here);
    
    if full||(i==i_here&&nargin==3)||(i==i_here&&j_here==0)
        J((i_sig-1)*6+1:i_sig*6,:) = J_here;
        i_sig                      = i_sig+1;
    end
    
    if Tr.VLinks(Tr.LinkIndex(i)).linktype == 'r'
        
        gi     = Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here = g_here*gi;
        J_here = dinamico_Adjoint(ginv(gi))*J_here;
        
        if full||(i==i_here&&nargin==3)
            J((i_sig-1)*6+1:i_sig*6,:) = J_here;
            i_sig                      = i_sig+1;
        end
        
        % bringing all quantities to the end of rigid link
        gf     = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here = g_here*gf;
        J_here = dinamico_Adjoint(ginv(gf))*J_here;
    end
    
    dof_start = dof_start+dof_here;
    
    for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link
        
        dof_here = Tr.CVTwists{i}(j+1).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);
        xi_star  = Tr.CVTwists{i}(j+1).xi_star;
        gi       = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        lpf      = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Xs       = Tr.CVTwists{i}(j+1).Xs;
        nip      = Tr.CVTwists{i}(j+1).nip;

        if Tr.Z_order==4
            B_Z1       = Tr.CVTwists{i}(j+1).B_Z1;
            B_Z2       = Tr.CVTwists{i}(j+1).B_Z2;
        else %order 2
            B_Z        = Tr.CVTwists{i}(j+1).B_Z;
        end

        
        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here                     = g_here*gi;
        J_here                     = dinamico_Adjoint(ginv(gi))*J_here;
        
        if full||(i==i_here&&nargin==3)||(i==i_here&&j==j_here)
            J((i_sig-1)*6+1:i_sig*6,:) = J_here;
            i_sig                      = i_sig+1;
        end
        
        Lscale         = lpf;
        g_here(1:3,4)  = g_here(1:3,4)/Lscale;
        J_here(4:6,:)  = J_here(4:6,:)/Lscale;
        
        
        for ii = 2:nip
            
            H    = Xs(ii)-Xs(ii-1);
            
            if Tr.Z_order==4
                
                xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1),2);
                xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
                xi_Z1here(1:3) = xi_Z1here(1:3)*Lscale; %scaling
                xi_Z2here(1:3) = xi_Z2here(1:3)*Lscale;
                
                B_Z1here  = B_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
                B_Z2here  = B_Z2(6*(ii-2)+1:6*(ii-1),:);
                if dof_here>0
                    xi_Z1here = B_Z1here*q_here+xi_Z1here;
                    xi_Z2here = B_Z2here*q_here+xi_Z2here;
                end
                ad_xi_Z1here = dinamico_adj(xi_Z1here);
                BGamma_here  = (H/2)*(B_Z1here+B_Z2here)+...
                               ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here); %choice of order from use    

                Gamma_here   = (H/2)*(xi_Z1here+xi_Z2here)+...
                                ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;
                      
            else % order 2
                xi_Zhere = xi_star(6*(ii-2)+1:6*(ii-1),4);
                xi_Zhere(1:3) = xi_Zhere(1:3)*Lscale; %scaling
                
                B_Zhere  = B_Z(6*(ii-2)+1:6*(ii-1),:);%note this step
                if dof_here>0
                    xi_Zhere = B_Zhere*q_here+xi_Zhere;
                end
                BGamma_here = H*B_Zhere;  

                
                Gamma_here  = H*xi_Zhere;

            end
            
            [gh,TGamma_here] = variable_expmap_gTg(Gamma_here);
            
            TBGamma_here = zeros(6,ndof);
            TBGamma_here(:,dof_start:dof_start+dof_here-1) = TGamma_here*BGamma_here;
            
            %updating g, Jacobian, Jacobian_dot and eta
            g_here = g_here*gh;
            J_here = dinamico_Adjoint(ginv(gh))*(J_here+TBGamma_here);
            
            J_heret                    = J_here;
            J_heret(4:6,:)             = J_heret(4:6,:)*Lscale; %back to SI units
            
            if full||(i==i_here&&nargin==3)||(i==i_here&&j==j_here)
                J((i_sig-1)*6+1:i_sig*6,:) = J_heret;
                i_sig                      = i_sig+1;
            end
            
        end
        
        g_here(1:3,4)  = g_here(1:3,4)*Lscale;
        J_here(4:6,:)  = J_here(4:6,:)*Lscale;
        
        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gf     = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
        g_here = g_here*gf;
        J_here = dinamico_Adjoint(ginv(gi))*J_here;
        
        dof_start = dof_start+dof_here;
    end
    g_Ltip((i-1)*4+1:i*4,:) = g_here;
    J_Ltip((i-1)*6+1:i*6,:) = J_here;
end
end

