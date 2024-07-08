%Function to display and save the actuation cables 
%Last modified by Anup Teejo Mathew (14.06.2021)
function CableDisplay(Tr,dc_fn,iL,div_start,div_end)

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

for i = 1:iL %upto iL is enough
    
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
        xi         = B_here*q_here+xi_star;
        g_joint    = variable_expmap_g(xi);
    end
    
    g_here = g_here*g_joint;

    if VLinks(LinkIndex(i)).linktype=='r'
        g_here  = g_here*VLinks(LinkIndex(i)).gi*VLinks(LinkIndex(i)).gf;
    end
    
    dof_start = dof_start+dof_here;
    %=============================================================================
    %Soft pieces
    for j = 1:(VLinks(LinkIndex(i)).npie)-1
        
        xi_starfn = VTwists(j+1).xi_starfn;
        gi      = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        Lscale  = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Type    = Tr.CVTwists{i}(j+1).Type;

        Bdof      = Tr.CVTwists{i}(j+1).Bdof;
        Bodr      = Tr.CVTwists{i}(j+1).Bodr;

        Bh     = Tr.CVTwists{i}(j+1).Bh;
        
        %updating g, Jacobian, Jacobian_dot and eta at X=0
        g_here                     = g_here*gi;
        dof_here  = VTwists(j+1).dof;
        q_here    = q(dof_start:dof_start+dof_here-1);

        n_l      = VLinks(LinkIndex(i)).n_l;
        if i==iL&&j>=div_start&&j<=div_end
            n_l=n_l*4; %more points
        end
        
        Xs       = linspace(0,1,n_l);
        H        = Xs(2)-Xs(1);
        Z = (1/2)*H;          % Zanna quadrature coefficient
        
        if j==div_start
            x_c=0;
        end
        if i==iL&&j==div_start&&j<=div_end
            dc_here   = dc_fn(x_c);
            PosC_here = g_here*[dc_here;1];
            XC        = [XC PosC_here(1)];
            YC        = [YC PosC_here(2)];
            ZC        = [ZC PosC_here(3)];
        end
        
        for ii = 1:n_l-1
            
            X   = Xs(ii);
            X_Z = X+Z;
            
            xi_Zhere  = xi_starfn(X_Z);
            xi_Zhere(1:3) = xi_Zhere(1:3)*Lscale;
            
            if ~isempty(q_here)                 
                if strcmp(Type,'FEM Like')
                    SubClass  = Tr.CVTwists{i}(j+1).SubClass;
                    xi_Zhere  = xi_Zhere+Bh(X_Z,Bdof,Bodr,SubClass)*q_here;
                elseif strcmp(Type,'Custom Independent')
                    xi_Zhere  = xi_Zhere+Bh(X_Z)*q_here;
                else
                    xi_Zhere  = xi_Zhere+Bh(X_Z,Bdof,Bodr)*q_here;
                end
            end
            Gamma_here = H*xi_Zhere;
            
            Gamma_here(4:6) = Gamma_here(4:6)*Lscale;
            gh         = variable_expmap_g(Gamma_here);
            g_here     = g_here*gh;
            
            if i==iL&&j>=div_start&&j<=div_end
                x_c       = x_c+H*Lscale;
                dc_here   = dc_fn(x_c);
                PosC_here = g_here*[dc_here;1];
                XC        = [XC PosC_here(1)];
                YC        = [YC PosC_here(2)];
                ZC        = [ZC PosC_here(3)];
            end

        end
        dof_start  = dof_start+dof_here;
    end
    g_Ltip((i-1)*4+1:i*4,:) = g_here; 
end

plot3(XC,YC,ZC,'LineWidth',2,'Color','m');

end