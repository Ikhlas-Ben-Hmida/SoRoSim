%Function that calculates the derivative of the jacobian with respect to
%time (Jd) at every significant point (24.05.2021)

function Jd = Jacobiandot(Tr,q,qd,i_here,j_here) %i_here is link index, j_here is division (0 for joint)

if isrow(q)
    q=q';
end
if isrow(qd)
    qd=qd';
end

N         = Tr.N;
ndof      = Tr.ndof;
g_ini     = Tr.g_ini; %initial configuration of all link wrt its previous link
iLpre     = Tr.iLpre;
g_Ltip    = repmat(eye(4),N,1);
Jd_Ltip   = repmat(zeros(6,ndof),N,1);

full = false;
if nargin==5
    if j_here==0
        nsig = 1;
    else
        nsig = Tr.CVTwists{i_here}(j_here+1).nip; %j_here>1 is allowed only for soft links
    end
elseif nargin==4
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
Jd        = zeros(6*nsig,ndof);
dof_start = 1; %starting dof of current piece
i_sig     = 1;


for i = 1:N
    
    if iLpre(i)>0
        g_here       = g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        Jd_here      = Ad_g_ini_inv*Jd_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        eta_here   = Ad_g_ini_inv*eta_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6);
    else
        g_here   = g_ini((i-1)*4+1:i*4,:);
        Jd_here  = zeros(6,ndof);
        eta_here   = zeros(6,1);
    end
    
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
        xi               = B_here*q_here+xi_star;
        xid              = B_here*qd_here;
        [g_joint,Tg,Tgd] = variable_expmap_gTgTgd_mex(xi,xid);

        TgB_here                                    = zeros(6,ndof);
        TgB_here(:,dof_start:dof_start+dof_here-1)  = Tg*B_here;
        TgBd_here                                   = zeros(6,ndof);
        TgBd_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*Tg*B_here+Tgd*B_here;
    end
    
    %updating g, Jacobian, Jacobian_dot and eta
    g_here         = g_here*g_joint;
    Ad_g_joint_inv = dinamico_Adjoint(ginv(g_joint));
    Jd_here        = Ad_g_joint_inv*(Jd_here+TgBd_here);
    eta_here       = Ad_g_joint_inv*(eta_here+TgB_here(:,dof_start:dof_start+dof_here-1)*qd_here);
    
    if full||(i==i_here&&nargin==4)||(i==i_here&&j_here==0)
        Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
        i_sig                       = i_sig+1;
    end
    
    if Tr.VLinks(Tr.LinkIndex(i)).linktype == 'r'
        
        gi        = Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here    = g_here*gi;
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        Jd_here   = Ad_gi_inv*Jd_here;
        eta_here  = Ad_gi_inv*eta_here;
        
        if full||(i==i_here&&nargin==4)
            Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
            i_sig                       = i_sig+1;
        end
        
        % bringing all quantities to the end of rigid link
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here    = g_here*gf;
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        Jd_here   = Ad_gf_inv*Jd_here;
        eta_here  = Ad_gf_inv*eta_here;
        
    end
    
    dof_start = dof_start+dof_here;
    
    for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link
        
        dof_here = Tr.CVTwists{i}(j+1).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);
        qd_here  = qd(dof_start:dof_start+dof_here-1);
        xi_star  = Tr.CVTwists{i}(j+1).xi_star;
        gi       = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        ld       = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Xs       = Tr.CVTwists{i}(j+1).Xs;
        nip      = Tr.CVTwists{i}(j+1).nip;
        
        
        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here   = g_here*gi;
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        Jd_here   = Ad_gi_inv*Jd_here;
        eta_here  = Ad_gi_inv*eta_here;
        
        if full||(i==i_here&&nargin==4)||(i==i_here&&j==j_here)
            Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_here;
            i_sig                       = i_sig+1;
        end 
        
        %scaling of quantities using the formula: Lscale m = 1 unit
        Lscale         = ld;
        g_here(1:3,4)  = g_here(1:3,4)/Lscale;
        Jd_here(4:6,:) = Jd_here(4:6,:)/Lscale;
        eta_here(4:6)  = eta_here(4:6)/Lscale;
        
        for ii = 2:nip
            
            H    = Xs(ii)-Xs(ii-1);
            
            if Tr.Z_order==4
                
                xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1),2); 
                xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
                xi_Z1here(1:3) = xi_Z1here(1:3)*Lscale; %scaling
                xi_Z2here(1:3) = xi_Z2here(1:3)*Lscale;
                
                %B is Phi for independent basis and Bq is Psi for dependent basis
                    
                B_Z1here  = Tr.CVTwists{i}(j+1).B_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
                B_Z2here  = Tr.CVTwists{i}(j+1).B_Z2(6*(ii-2)+1:6*(ii-1),:);

                if dof_here>0
                    xi_Z1here = B_Z1here*q_here+xi_Z1here;
                    xi_Z2here = B_Z2here*q_here+xi_Z2here;

                    xid_Z1here  = B_Z1here*qd_here;

                    ad_xi_Z1here = dinamico_adj(xi_Z1here);

                    BGamma_here  = (H/2)*(B_Z1here+B_Z2here)+...
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
            Jd_here    = Ad_gh_inv*(Jd_here+TBGammad_here); %full
            eta_here = Ad_gh_inv*(eta_here+TBGamma_here(:,dof_start:dof_start+dof_here-1)*qd_here);

            
            Jd_heret                    = Jd_here;
            Jd_heret(4:6,:)             = Jd_heret(4:6,:)*Lscale;
            
            if full||(i==i_here&&nargin==4)||(i==i_here&&j==j_here)
                Jd((i_sig-1)*6+1:i_sig*6,:) = Jd_heret;
                i_sig                       = i_sig+1;
            end
            
        end
        
        g_here(1:3,4)  = g_here(1:3,4)*Lscale;
        Jd_here(4:6,:) = Jd_here(4:6,:)*Lscale;
        eta_here(4:6)  = eta_here(4:6)*Lscale;
        
        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
        g_here    = g_here*gf;
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        Jd_here   = Ad_gf_inv*Jd_here;
        eta_here  = Ad_gf_inv*eta_here;
        
        dof_start  = dof_start+dof_here;
    end
    g_Ltip((i-1)*4+1:i*4,:)   = g_here;
    Jd_Ltip((i-1)*6+1:i*6,:)  = Jd_here;
    eta_Ltip((i-1)*6+1:i*6,:) = eta_here;
end
end

