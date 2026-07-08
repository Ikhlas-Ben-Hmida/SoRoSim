function plotq(Linkage,q)

if nargin==1
    q=zeros(Linkage.ndof,1);
end

if isrow(q)
    q=q';
end

PlotParameters = Linkage.PlotParameters;

if PlotParameters.ClosePrevious
    close all
end

%Plot options

fh=figure(1);
fh.Units='normalized';
fh.OuterPosition=[0 0 1 1];

set(gca,'CameraPosition',PlotParameters.CameraPosition,...
    'CameraTarget',PlotParameters.CameraTarget,...
    'CameraUpVector',PlotParameters.CameraUpVector,...
    'FontSize',18)

if PlotParameters.Light
    camlight(PlotParameters.Az_light,PlotParameters.El_light)
end
% view(0,90)
axis equal
grid on
hold on
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

% Set all text elements to use LaTeX interpreter
set(get(gca, 'Title'), 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'Interpreter', 'latex');
set(get(gca, 'YLabel'), 'Interpreter', 'latex');
set(get(gca, 'ZLabel'), 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

set(gca,'FontSize',12)

axis ([PlotParameters.XLim PlotParameters.YLim PlotParameters.ZLim]);

%Forward Kinematics: (Product of exponentials)
N         = Linkage.N;
dof_start = 1;
g_ini     = Linkage.g_ini;
g_Ltip    = repmat(eye(4),N,1);
iLpre     = Linkage.iLpre;

for i=1:N
    
    if iLpre(i)>0
        g_here=g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
    else
        g_here=g_ini((i-1)*4+1:i*4,:);
    end
    
    %Rigid link or joint
    dof_here   = Linkage.CVRods{i}(1).dof;
    q_here     = q(dof_start:dof_start+dof_here-1);
    Phi_here   = Linkage.CVRods{i}(1).Phi;
    xi_star    = Linkage.CVRods{i}(1).xi_star;
    
    if dof_here == 0 %fixed joint (N)
        g_joint = eye(4);
    else
        xi         = Phi_here*q_here+xi_star;
        g_joint    = variable_expmap_g(xi);
    end
    g_here      = g_here*g_joint;
    
    n_r   = Linkage.VLinks(Linkage.LinkIndex(i)).n_r;
    if strcmp(Linkage.VLinks(Linkage.LinkIndex(i)).CS,'R')
        n_r=5;
    end

    n_l   = Linkage.VLinks(Linkage.LinkIndex(i)).n_l;
    color = Linkage.VLinks(Linkage.LinkIndex(i)).color;
    alpha = Linkage.VLinks(Linkage.LinkIndex(i)).alpha;

    if Linkage.VLinks(Linkage.LinkIndex(i)).L>0
    if Linkage.VLinks(Linkage.LinkIndex(i)).linktype=='r'
        
        L       = Linkage.VLinks(Linkage.LinkIndex(i)).L;
        gi      = Linkage.VLinks(Linkage.LinkIndex(i)).gi;
        g_here  = g_here*gi;
        
        if ~Linkage.VLinks(Linkage.LinkIndex(i)).CPF
            Xr      = linspace(0,L,n_l);
            g_hereR = g_here*[eye(3) [-Linkage.VLinks(Linkage.LinkIndex(i)).cx;0;0];0 0 0 1];
            dx      = Xr(2)-Xr(1);

            Xpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
            Ypatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
            Zpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
            i_patch = 1;

            [y,z] = computeBoundaryYZ(Linkage.VLinks(Linkage.LinkIndex(i)),0);
            pos  = [zeros(1,n_r);y;z;ones(1,n_r)]; %homogeneous positions in local frame 4xn_r

            pos_here = g_hereR*pos;
            x_here   = pos_here(1,:);
            y_here   = pos_here(2,:);
            z_here   = pos_here(3,:);

            Xpatch(:,i_patch) = x_here';
            Ypatch(:,i_patch) = y_here';
            Zpatch(:,i_patch) = z_here';
            i_patch           = i_patch+1;

            x_pre    = x_here;
            y_pre    = y_here;
            z_pre    = z_here;


            for ii=2:n_l

                [y,z] = computeBoundaryYZ(Linkage.VLinks(Linkage.LinkIndex(i)),Xr(ii)/L);
                pos  = [zeros(1,n_r);y;z;ones(1,n_r)]; %homogeneous positions in local frame 4xn_r

                g_hereR  = g_hereR*[eye(3) [dx;0;0];0 0 0 1];
                pos_here = g_hereR*pos;
                x_here   = pos_here(1,:);
                y_here   = pos_here(2,:);
                z_here   = pos_here(3,:);

                %Plotting rigid link
                for jj=1:n_r-1
                    Xpatch(1:5,i_patch)   = [x_pre(jj) x_here(jj) x_here(jj+1) x_pre(jj+1) x_pre(jj)]';
                    Xpatch(6:end,i_patch) = x_pre(jj)*ones(n_r-5,1);
                    Ypatch(1:5,i_patch)   = [y_pre(jj) y_here(jj) y_here(jj+1) y_pre(jj+1) y_pre(jj)]';
                    Ypatch(6:end,i_patch) = y_pre(jj)*ones(n_r-5,1);
                    Zpatch(1:5,i_patch)   = [z_pre(jj) z_here(jj) z_here(jj+1) z_pre(jj+1) z_pre(jj)]';
                    Zpatch(6:end,i_patch) = z_pre(jj)*ones(n_r-5,1);
                    i_patch = i_patch+1;
                end

                x_pre    = x_here;
                y_pre    = y_here;
                z_pre    = z_here;

            end

            Xpatch(:,i_patch) = x_here';
            Ypatch(:,i_patch) = y_here';
            Zpatch(:,i_patch) = z_here';

            patch(Xpatch,Ypatch,Zpatch,color,'EdgeColor','none','FaceAlpha',alpha);
        else
            CustomShapePlot(g_here);
        end
        gf     = Linkage.VLinks(Linkage.LinkIndex(i)).gf;
        g_here = g_here*gf;
    end
    end

    if ~Linkage.OneBasis
        dof_start = dof_start+dof_here;
    end
    
    %=============================================================================
    %Soft link pieces
    for j=1:(Linkage.VLinks(Linkage.LinkIndex(i)).npie)-1
        
        dof_here   = Linkage.CVRods{i}(j+1).dof;
        Type       = Linkage.CVRods{i}(j+1).Type;
        q_here     = q(dof_start:dof_start+dof_here-1);
        xi_starfn  = Linkage.CVRods{i}(j+1).xi_starfn;
        gi         = Linkage.VLinks(Linkage.LinkIndex(i)).gi{j};
        Phi_dof    = Linkage.CVRods{i}(j+1).Phi_dof;
        Phi_odr    = Linkage.CVRods{i}(j+1).Phi_odr;
        Phi_h      = Linkage.CVRods{i}(j+1).Phi_h;
        ld         = Linkage.VLinks(Linkage.LinkIndex(i)).ld{j};

        g_here = g_here*gi;

        Xs = linspace(0,1,n_l);
        H  = Xs(2)-Xs(1);
        
        Z = (1/2)*H;          % Zanna quadrature coefficient
        
        Xpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
        Ypatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
        Zpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
        i_patch = 1;
        
        %cross sectional shape Circular, Rectangular, and Ellipsoidal
        [y,z] = computeBoundaryYZ(Linkage.VLinks(Linkage.LinkIndex(i)),0,j);
        pos  = [zeros(1,n_r);y;z;ones(1,n_r)]; %homogeneous positions in local frame 4xn_r
        
        pos_here = g_here*pos;
        x_here   = pos_here(1,:);
        y_here   = pos_here(2,:);
        z_here   = pos_here(3,:);
        
        Xpatch(:,i_patch) = x_here';
        Ypatch(:,i_patch) = y_here';
        Zpatch(:,i_patch) = z_here';
        i_patch           = i_patch+1;
        
        x_pre = x_here;
        y_pre = y_here;
        z_pre = z_here;
        
        for ii=1:n_l-1
            %cross sectional shape Circular, Rectangular, and Ellipsoidal
            [y,z] = computeBoundaryYZ(Linkage.VLinks(Linkage.LinkIndex(i)),Xs(ii+1),j);
            pos  = [zeros(1,n_r);y;z;ones(1,n_r)]; %homogeneous positions in local frame 4xn_r
            
            X = Xs(ii);
            X_Z = X+Z;

            xi_Zhere  = xi_starfn(X_Z);
            Phi_Scale = diag([1/ld 1/ld 1/ld 1 1 1]);
            
            if ~isempty(q_here)
                if strcmp(Type,'FEM Like')
                    SubClass  = Linkage.CVRods{i}(j+1).SubClass;
                    xi_Zhere  = xi_Zhere+Phi_Scale*Phi_h(X_Z,Phi_dof,Phi_odr,SubClass)*q_here;
                elseif strcmp(Type,'Custom Independent')
                    xi_Zhere  = xi_Zhere+Phi_Scale*Phi_h(X_Z)*q_here;
                else
                    xi_Zhere  = xi_Zhere+Phi_Scale*Phi_h(X_Z,Phi_dof,Phi_odr)*q_here;
                end
            end
            
            Gamma_here = H*ld*xi_Zhere;

            gh         = variable_expmap_g(Gamma_here);
            g_here     = g_here*gh;

            pos_here = g_here*pos;
            x_here     = pos_here(1,:);
            y_here     = pos_here(2,:);
            z_here     = pos_here(3,:);
            
            %Plotting soft link pieces
            for jj=1:n_r-1
                Xpatch(1:5,i_patch)   = [x_pre(jj) x_here(jj) x_here(jj+1) x_pre(jj+1) x_pre(jj)]';
                Xpatch(6:end,i_patch) = x_pre(jj)*ones(n_r-5,1);
                Ypatch(1:5,i_patch)   = [y_pre(jj) y_here(jj) y_here(jj+1) y_pre(jj+1) y_pre(jj)]';
                Ypatch(6:end,i_patch) = y_pre(jj)*ones(n_r-5,1);
                Zpatch(1:5,i_patch)   = [z_pre(jj) z_here(jj) z_here(jj+1) z_pre(jj+1) z_pre(jj)]';
                Zpatch(6:end,i_patch) = z_pre(jj)*ones(n_r-5,1);
                i_patch = i_patch+1;
            end
            x_pre    = x_here;
            y_pre    = y_here;
            z_pre    = z_here;
        end
        
        Xpatch(:,i_patch) = x_here';
        Ypatch(:,i_patch) = y_here';
        Zpatch(:,i_patch) = z_here';

        patch(Xpatch,Ypatch,Zpatch,color,'EdgeColor','none','FaceAlpha',alpha);
        
        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gf     = Linkage.VLinks(Linkage.LinkIndex(i)).gf{j};
        g_here = g_here*gf;

        if ~Linkage.OneBasis
            dof_start = dof_start+dof_here;
        end
    end
    g_Ltip((i-1)*4+1:i*4,:) = g_here;
end
% drawnow
end