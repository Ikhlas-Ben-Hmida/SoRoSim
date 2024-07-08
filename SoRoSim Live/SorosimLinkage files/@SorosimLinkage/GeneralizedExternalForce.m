%Computes the generalized external force
%Last modified by Anup Teejo Mathew - 24/05/2021
function F = GeneralizedExternalForce(Tr,q,t,qd)

ndof    = Tr.ndof;
if nargin==1
    t=0;
    qd = zeros(ndof,1);
elseif nargin==2
    qd = zeros(ndof,1);
end

if isrow(q)
    q=q';
end

N       = Tr.N;
nsig    = Tr.nsig;
G       = Tr.G;
g_ini   = Tr.g_ini; %initial configuration of all link wrt its previous link
iLpre   = Tr.iLpre;
g_Ltip  = repmat(eye(4),N,1);
J_Ltip  = repmat(zeros(6,ndof),N,1);
g       = zeros(4*nsig,4);
J       = zeros(6*nsig,ndof);
eta_Ltip   = zeros(N*6,1); %total velocity J*qd
F   = zeros(ndof,1);%external force 
eta = zeros(6*nsig,1);

i_sig        = 1;
i_sig_nj     = 1;
dof_start    = 1; %starting dof of current piece
for i=1:N

    if iLpre(i)>0
        g_here       = g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        J_here       = Ad_g_ini_inv*J_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        eta_here     = Ad_g_ini_inv*eta_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6);

    else
        g_here   = g_ini((i-1)*4+1:i*4,:);
        J_here   = zeros(6,ndof);
        eta_here = zeros(6,1);
    end

    %Joint
    dof_here = Tr.CVTwists{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    qd_here  = q(dof_start:dof_start+dof_here-1);
    B_here   = Tr.CVTwists{i}(1).B;
    xi_star  = Tr.CVTwists{i}(1).xi_star;

    if dof_here==0 %fixed joint (N)
        g_joint   = eye(4);
        TgB_here  = zeros(6,ndof);
    else
        xi           = B_here*q_here+xi_star;
        [g_joint,Tg] = variable_expmap_gTg(xi);

        TgB_here                                    = zeros(6,ndof);
        TgB_here(:,dof_start:dof_start+dof_here-1)  = Tg*B_here;
    end

    %updating g, Jacobian, Jacobian_dot and eta
    g_here         = g_here*g_joint;
    Ad_g_joint_inv = dinamico_Adjoint(ginv(g_joint));
    J_here         = Ad_g_joint_inv*(J_here+TgB_here);
    eta_here       = Ad_g_joint_inv*(eta_here+TgB_here(:,dof_start:dof_start+dof_here-1)*qd_here);

    g((i_sig-1)*4+1:i_sig*4,:) = g_here;
    J((i_sig-1)*6+1:i_sig*6,:) = J_here;
    eta((i_sig-1)*6+1:i_sig*6) = eta_here;

    i_sig = i_sig+1;
    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'

        gi     = Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here = g_here*gi;
        J_here = dinamico_Adjoint(ginv(gi))*J_here;
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        J_here    = Ad_gi_inv*J_here;
        eta_here  = Ad_gi_inv*eta_here;

        g((i_sig-1)*4+1:i_sig*4,:) = g_here;
        J((i_sig-1)*6+1:i_sig*6,:) = J_here;
        eta((i_sig-1)*6+1:i_sig*6)  = eta_here;

        M_here=Tr.VLinks(Tr.LinkIndex(i)).M;
        if Tr.Gravity
            F = F+J_here'*M_here*dinamico_Adjoint(ginv(g_here))*G;
        end

        % bringing all quantities to the end of rigid link
        gf         = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here     = g_here*gf;
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        J_here    = Ad_gf_inv*J_here;
        eta_here  = Ad_gf_inv*eta_here;

        i_sig    = i_sig+1;
        i_sig_nj = i_sig_nj+1;

    end

    dof_start=dof_start+dof_here;

    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link

        dof_here   = Tr.CVTwists{i}(j+1).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        qd_here    = qd(dof_start:dof_start+dof_here-1);

        xi_star = Tr.CVTwists{i}(j+1).xi_star;
        gi      = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        lpf     = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Ms      = Tr.CVTwists{i}(j+1).Ms;
        Xs      = Tr.CVTwists{i}(j+1).Xs;
        Ws      = Tr.CVTwists{i}(j+1).Ws;
        nip     = Tr.CVTwists{i}(j+1).nip;

        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here    = g_here*gi;
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        J_here    = Ad_gi_inv*J_here;
        eta_here  = Ad_gi_inv*eta_here;

        %scaling of quantities
        Lscale         = lpf;
        g_here(1:3,4)  = g_here(1:3,4)/Lscale;
        J_here(4:6,:)  = J_here(4:6,:)/Lscale;
        G = G/Lscale;

        g((i_sig-1)*4+1:i_sig*4,:)  = g_here;
        J((i_sig-1)*6+1:i_sig*6,:)  = J_here;
        eta((i_sig-1)*6+1:i_sig*6)  = eta_here;

        ii = 1;
        if Ws(ii)>0
            W_here  = Ws(ii);
            Ms_here = Ms(6*(ii-1)+1:6*ii,:);

            %scaling
            Ms_here(1:3,:) = Ms_here(1:3,:)/Lscale;
            Ms_here(4:6,:) = Ms_here(4:6,:)*Lscale;
            if Tr.Gravity
                F = F+W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %rescale
            end
        end

        i_sig                       = i_sig+1;
        i_sig_nj                    = i_sig_nj+1;

        for ii=2:nip

            H = Xs(ii)-Xs(ii-1);

            if Tr.Z_order==4

                xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1),2); 
                xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
                xi_Z1here(1:3) = xi_Z1here(1:3)*Lscale; %scaling
                xi_Z2here(1:3) = xi_Z2here(1:3)*Lscale;

                B_Z1here  = Tr.CVTwists{i}(j+1).B_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
                B_Z2here  = Tr.CVTwists{i}(j+1).B_Z2(6*(ii-2)+1:6*(ii-1),:);

                if dof_here>0

                    xi_Z1here = B_Z1here*q_here+xi_Z1here;
                    xi_Z2here = B_Z2here*q_here+xi_Z2here;

                    ad_xi_Z1here = dinamico_adj(xi_Z1here);

                    BGamma_here  = (H/2)*(B_Z1here+B_Z2here)+...
                                       ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);


                else

                    ad_xi_Z1here = dinamico_adj(xi_Z1here);
                    BGamma_here  = (H/2)*(B_Z1here+B_Z2here)+...
                                       ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);

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
                else
                    BGamma_here = H*B_Zhere;
                end

                Gamma_here  = H*xi_Zhere;

            end

            [gh,TGamma_here] = variable_expmap_gTg(Gamma_here);

            TBGamma_here                                   = zeros(6,ndof);
            TBGamma_here(:,dof_start:dof_start+dof_here-1) = TGamma_here*BGamma_here;

            %updating g, Jacobian, Jacobian_dot and eta
            g_here     = g_here*gh;
            Ad_gh_inv  = dinamico_Adjoint(ginv(gh));
            J_here     = Ad_gh_inv*(J_here+TBGamma_here);
            eta_here   = Ad_gh_inv*(eta_here+TBGamma_here(:,dof_start:dof_start+dof_here-1)*qd_here);


            g((i_sig-1)*4+1:i_sig*4,:)  = g_here;
            J((i_sig-1)*6+1:i_sig*6,:)  = J_here;
            eta((i_sig-1)*6+1:i_sig*6)  = eta_here;

            %integrals evaluation
            if Ws(ii)>0
                W_here  = Ws(ii);
                Ms_here = Ms(6*(ii-1)+1:6*ii,:);

                %scaling
                Ms_here(1:3,:) = Ms_here(1:3,:)/Lscale;
                Ms_here(4:6,:) = Ms_here(4:6,:)*Lscale;
                if Tr.Gravity
                    F = F+W_here*J_here'*Ms_here*dinamico_Adjoint(ginv(g_here))*G*Lscale^2; %rescale
                end
 

            end
            i_sig    = i_sig+1;
            i_sig_nj = i_sig_nj+1;
        end
        %scaling back quantities
        g_here(1:3,4)  = g_here(1:3,4)*Lscale;
        J_here(4:6,:)  = J_here(4:6,:)*Lscale;
        eta_here(4:6)  = eta_here(4:6)*Lscale;

        G              = G*Lscale;

        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
        g_here    = g_here*gf;
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        J_here    = Ad_gf_inv*J_here;
        eta_here  = Ad_gf_inv*eta_here;

        dof_start = dof_start+dof_here;
    end
    g_Ltip((i-1)*4+1:i*4,:)  = g_here;
    J_Ltip((i-1)*6+1:i*6,:)  = J_here;
    eta_Ltip((i-1)*6+1:i*6,:)  = eta_here;

end
%% Point Force
if Tr.PointForce
    F = F+ComputePointForce(Tr,J,g,t);
end

%scale g, J, and eta
if Tr.CEFP||Tr.CAP

    i_sig    = 1;
    g_act    = g;
    J_act    = J;
    eta_act  = eta;

    for i=1:N
        i_sig = i_sig+1;
        if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
            i_sig = i_sig+1;
        end
        for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
            Lscale = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
            for ii=1:Tr.CVTwists{i}(j+1).nip
                g_here        = g((i_sig-1)*4+1:i_sig*4,:);
                g_here(1:3,4) = g_here(1:3,4)*Lscale;
                g_act((i_sig-1)*4+1:i_sig*4,:) = g_here;

                J_here        = J((i_sig-1)*6+1:i_sig*6,:);
                J_here(4:6,:) = J_here(4:6,:)*Lscale;
                J_act((i_sig-1)*6+1:i_sig*6,:) = J_here;

                eta_here      = eta((i_sig-1)*6+1:i_sig*6);
                eta_here(4:6) = eta_here(4:6,:)*Lscale;
                eta_act((i_sig-1)*6+1:i_sig*6) = eta_here;

                i_sig = i_sig+1;               
            end
        end
    end

    if Tr.CEFP
        Fext = CustomExtForce(Tr,q,g_act,J_act,t,qd,eta_act,[]);%Jd not included
        FextP = CustomExtPointForce(Tr,q,g_act,J_act,t,qd,eta_act,[]);
        F_custom = CustomExtForce_Qspace(Tr,J,Fext,FextP);
        F = F+F_custom;
    end
    if Tr.CAP
        Fact = CustomActuation(Tr,q,g_act,J_act,t,qd,eta_act,[]);%Jd not included
        TauC = CustomActuation_Qspace(Tr,Fact,q,t);
        F    = F+TauC; %added with F
    end

end

end

