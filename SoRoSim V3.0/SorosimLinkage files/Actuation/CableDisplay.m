%Function to display and save the actuation cables 
%Last modified by Anup Teejo Mathew (14.06.2021)
function CableDisplay(y_fn,z_fn,i_sact,iL,div_start,div_end,Tr,done)
% 
% close all
% S.plotq0;

%Forward kinematics
q         = zeros(Tr.ndof,1);
N         = Tr.N;
dof_start = 1;
g_ini     = Tr.g_ini;
g_Ltip    = repmat(eye(4),N,1);
iLpre     = Tr.iLpre;
VLinks    = Tr.VLinks;
CVTwists  = Tr.CVTwists;
LinkIndex = Tr.LinkIndex;

XC        = [];
YC        = [];
ZC        = [];

for i = 1:iL
    
    if iLpre(i)>0
        g_here = g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
    else
        g_here = g_ini((i-1)*4+1:i*4,:);
    end
    
    VTwists = CVTwists{i};
    
    %Joints
    xi_star  = VTwists(1).xi_star;
    B_here   = VTwists(1).B;
    dof_here = VTwists(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    
    if dof_here == 0 %fixed joint (N)
        g_joint = eye(4);
    else
        if VLinks(LinkIndex(i)).jointtype=='U' %special case for universal joint. Considered as 2 revolute joints
            % first revolute joint
            xi         = B_here(:,1)*q_here(1)+xi_star;
            g_joint    = joint_expmap(xi);
            g_here     = g_here*g_joint;

            % second revolute joint
            xi         = B_here(:,2)*q_here(2)+xi_star;
            g_joint    = joint_expmap(xi);
        else
            xi         = B_here*q_here+xi_star;
            g_joint    = joint_expmap(xi);
        end
    end
    
    g_here = g_here*g_joint;

    if VLinks(LinkIndex(i)).linktype=='r'
        g_here  = g_here*[eye(3) [VLinks(LinkIndex(i)).L;0;0];0 0 0 1];
    end
    
    dof_start = dof_start+dof_here;
    %=============================================================================
    %Soft pieces
    for j = 1:(VLinks(LinkIndex(i)).npie)-1
        
        xi_starfn = VTwists(j+1).xi_starfn;
        gi        = VLinks(LinkIndex(i)).gi{j};
        Bdof      = VTwists(j+1).Bdof;
        Bodr      = VTwists(j+1).Bodr;
        ld        = VLinks(LinkIndex(i)).lp{j};
        g_here    = g_here*gi;
        dof_here  = VTwists(j+1).dof;
        q_here    = q(dof_start:dof_start+dof_here-1);

        n_l      = VLinks(LinkIndex(i)).n_l;
        if i==iL&&j>=div_start&&j<=div_end
            n_l=n_l*4;
        end
        
        if j==div_start
            l_pre=0;
        end
        
        Xs       = linspace(0,ld,n_l);
        H        = Xs(2)-Xs(1);
        Z1       = (1/2-sqrt(3)/6)*H;      %Zanna quadrature coefficient
        Z2       = (1/2+sqrt(3)/6)*H;      %Zanna quadrature coefficient
        B_Z1here = zeros(6,dof_here);
        B_Z2here = zeros(6,dof_here);
        
        if i==iL&&j>=div_start&&j<=div_end
            x_here    = l_pre;
            y_here    = y_fn(x_here);
            z_here    = z_fn(x_here);

            PosC_here = g_here*[0;y_here;z_here;1];
            XC        = [XC PosC_here(1)];
            YC        = [YC PosC_here(2)];
            ZC        = [ZC PosC_here(3)];
        end
        
        for ii = 1:n_l-1
            
            X    = Xs(ii);
            X_Z1 = X+Z1;
            X_Z2 = X+Z2;
            
            for jj=1:6
                for k = 1:Bdof(jj)*Bodr(jj)+Bdof(jj)
                    kk              = sum(Bdof(1:jj-1).*Bodr(1:jj-1))+sum(Bdof(1:jj-1))+k;
                    B_Z1here(jj,kk) = X_Z1^(k-1);
                    B_Z2here(jj,kk) = X_Z2^(k-1);
                end
            end

            if ~isempty(q_here)
                xi_Z1here  = B_Z1here*q_here+xi_starfn(X_Z1);
                xi_Z2here  = B_Z2here*q_here+xi_starfn(X_Z2);
            else
                xi_Z1here  = xi_starfn(X_Z1);
                xi_Z2here  = xi_starfn(X_Z2);
            end
            
            Gamma_here = (H/2)*(xi_Z1here+xi_Z2here)+...
                ((sqrt(3)*H^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
            
            k_here     = Gamma_here(1:3);
            theta_here = norm(k_here);
            gh         = variable_expmap(theta_here,Gamma_here);
            g_here     = g_here*gh;
            
            if i==iL&&j>=div_start&&j<=div_end
                x_here    = l_pre+Xs(ii+1);
                y_here    = y_fn(x_here);
                z_here    = z_fn(x_here);
                PosC_here = g_here*[0;y_here;z_here;1];
                XC        = [XC PosC_here(1)];
                YC        = [YC PosC_here(2)];
                ZC        = [ZC PosC_here(3)];
            end

        end
        dof_start  = dof_start+dof_here;
        if j>=div_start&&div_start>0
            l_pre = l_pre+ld;
        end
    end
    g_Ltip((i-1)*4+1:i*4,:) = g_here; 
end

plot3(XC,YC,ZC,'LineWidth',2,'Color','m');

if done
    
    if exist('CablePoints.mat', 'file')
        load('CablePoints.mat','Cy_fn','Cz_fn');
        Cy_fn{i_sact,iL} = y_fn;
        Cz_fn{i_sact,iL} = z_fn;
        save('CablePoints.mat','Cy_fn','Cz_fn');
    else
        Cy_fn = cell(Tr.n_sact,Tr.N); 
        Cz_fn = cell(Tr.n_sact,Tr.N);
        Cy_fn{i_sact,iL} = y_fn;
        Cz_fn{i_sact,iL} = z_fn;
        save('CablePoints.mat','Cy_fn','Cz_fn');
    end
    
end

end