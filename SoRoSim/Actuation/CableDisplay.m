%Function to display and save the actuation cables 
%Last modified by Anup Teejo Mathew (14.06.2021)
function CableDisplay(y_sym,z_sym,iL,Sdiv,Ediv,S,done)
% 
% close all
% S.plotq0;

%Forward kinematics
q         = zeros(S.ndof,1);
N         = S.N;
f         = 1;
dof_start = 1;
g_here    = S.g_ini;

XC        = [];
YC        = [];
ZC        = [];

for i = 1:iL
        
    %Joints
    xi_star  = S.Vtwists(f).xi_star;
    B_here   = S.Vtwists(f).B;
    dof_here = S.Vtwists(f).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    
    if dof_here == 0 %fixed joint (N)
        g_joint = eye(4);
    else
         if S.VLinks(i).jointtype=='U' %special case for universal joint. Considered as 2 revolute joints
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

    if S.VLinks(i).linktype=='r'
        g_here  = g_here*[eye(3) [S.VLinks(i).L;0;0];0 0 0 1];
    end
    
    f=f+1;
    dof_start = dof_start+dof_here;
    %=============================================================================
    %Soft pieces
    for j = 1:(S.VLinks(i).npie)-1
        
        xi_starfn = S.Vtwists(f).xi_starfn;
        Bdof      = S.Vtwists(f).Bdof;
        Bodr      = S.Vtwists(f).Bodr;
        lpf       = S.VLinks(i).lp{j};
        dof_here  = S.Vtwists(f).dof;
        q_here    = q(dof_start:dof_start+dof_here-1);

        n_l      = S.VLinks(i).n_l;
        if i==iL&&j>=Sdiv&&j<=Ediv
            n_l=n_l*4;
        end
        
        if j==Sdiv
            l_pre=0;
        end
        
        xS       = linspace(0,lpf,n_l);
        H        = xS(2)-xS(1);
        Z1       = (1/2-sqrt(3)/6)*H;      %Zanna quadrature coefficient
        Z2       = (1/2+sqrt(3)/6)*H;      %Zanna quadrature coefficient
        B_Z1here = zeros(6,dof_here);
        B_Z2here = zeros(6,dof_here);
        
        if i==iL&&j>=Sdiv&&j<=Ediv
            x_here    = l_pre;
            y_here    = double(subs(y_sym,x_here));
            z_here    = double(subs(z_sym,x_here));

            PosC_here = g_here*[0;y_here;z_here;1];
            XC        = [XC PosC_here(1)];
            YC        = [YC PosC_here(2)];
            ZC        = [ZC PosC_here(3)];
        end
        
        for ii = 1:n_l-1
            
            x    = xS(ii);
            x_Z1 = x+Z1;
            x_Z2 = x+Z2;
            
            for jj=1:6
                for k = 1:Bdof(jj)*Bodr(jj)+Bdof(jj)
                    kk              = sum(Bdof(1:jj-1).*Bodr(1:jj-1))+sum(Bdof(1:jj-1))+k;
                    B_Z1here(jj,kk) = x_Z1^(k-1);
                    B_Z2here(jj,kk) = x_Z2^(k-1);
                end
            end
            
            if ~isempty(q_here)
                xi_Z1here  = B_Z1here*q_here+xi_starfn(x_Z1);
                xi_Z2here  = B_Z2here*q_here+xi_starfn(x_Z2);
            else
                xi_Z1here  = xi_starfn(x_Z1);
                xi_Z2here  = xi_starfn(x_Z2);
            end
            
            Gamma_here = (H/2)*(xi_Z1here+xi_Z2here)+...
                ((sqrt(3)*H^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
            
            k_here     = Gamma_here(1:3);
            theta_here = norm(k_here);
            gh         = variable_expmap(theta_here,Gamma_here);
            g_here     = g_here*gh;
            
            if i==iL&&j>=Sdiv&&j<=Ediv
                x_here    = l_pre+xS(ii+1);
                y_here    = double(subs(y_sym,x_here));
                z_here    = double(subs(z_sym,x_here));
                PosC_here = g_here*[0;y_here;z_here;1];
                XC        = [XC PosC_here(1)];
                YC        = [YC PosC_here(2)];
                ZC        = [ZC PosC_here(3)];
            end

        end
        dof_start  = dof_start+dof_here;
        f=f+1;
        if j>=Sdiv
            l_pre = l_pre+lpf;
        end
    end
end

plot3(XC,YC,ZC,'LineWidth',2,'Color','m');

if done
    
    if exist('CablePoints.mat', 'file')
        
      load('CablePoints.mat','XC_save','YC_save','ZC_save');
      XC_save=[XC_save XC NaN];
      YC_save=[YC_save YC NaN];
      ZC_save=[ZC_save ZC NaN];
      
      save('CablePoints.mat','XC_save','YC_save','ZC_save');
    else
       XC_save = [XC NaN]; %for the next condition
       YC_save = [YC NaN];
       ZC_save = [ZC NaN];
      save('CablePoints.mat','XC_save','YC_save','ZC_save');
    end
    
end

end