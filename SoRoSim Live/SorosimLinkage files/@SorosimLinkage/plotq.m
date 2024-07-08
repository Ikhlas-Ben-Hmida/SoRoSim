function plotq(Tr,q)

if nargin==1
    q=zeros(Tr.ndof,1);
end

if isrow(q)
    q=q';
end

PlottingParameters = Tr.PlotParameters;

if PlottingParameters.ClosePrevious
    close all
end

%Plot options

fh=figure(1);
fh.Units='normalized';
fh.OuterPosition=[0 0 1 1];

set(gca,'CameraPosition',PlottingParameters.CameraPosition,...
    'CameraTarget',PlottingParameters.CameraTarget,...
    'CameraUpVector',PlottingParameters.CameraUpVector,...
    'FontSize',18)

if PlottingParameters.Light
    camlight(PlottingParameters.Az_light,PlottingParameters.El_light)
end
% view(0,90)
axis equal
grid on
hold on
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

axis ([PlottingParameters.X_lim PlottingParameters.Y_lim PlottingParameters.Z_lim]);

%Forward Kinematics: (Product of exponentials)
N         = Tr.N;
dof_start = 1;
g_ini     = Tr.g_ini;
g_Ltip    = repmat(eye(4),N,1);
iLpre     = Tr.iLpre;

for i=1:N
    
    if iLpre(i)>0
        g_here=g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
    else
        g_here=g_ini((i-1)*4+1:i*4,:);
    end
    
    %Rigid link or joint
    dof_here   = Tr.CVTwists{i}(1).dof;
    q_here     = q(dof_start:dof_start+dof_here-1);
    B_here     = Tr.CVTwists{i}(1).B;
    xi_star    = Tr.CVTwists{i}(1).xi_star;
    
    if dof_here == 0 %fixed joint (N)
        g_joint = eye(4);
    else
        xi         = B_here*q_here+xi_star;
        g_joint    = variable_expmap_g(xi);
    end
    g_here      = g_here*g_joint;
    
    n_r   = Tr.VLinks(Tr.LinkIndex(i)).n_r;
    if Tr.VLinks(Tr.LinkIndex(i)).CS=='R'
        n_r=5;
    end
    n_l   = Tr.VLinks(Tr.LinkIndex(i)).n_l;
    color = Tr.VLinks(Tr.LinkIndex(i)).color;
    if Tr.VLinks(Tr.LinkIndex(i)).L>0
    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        
        L       = Tr.VLinks(Tr.LinkIndex(i)).L;
        gi      = Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here  = g_here*gi;
        
        if ~Tr.VLinks(Tr.LinkIndex(i)).CPF
            Xr      = linspace(0,L,n_l);
            g_hereR = g_here*[eye(3) [-Tr.VLinks(Tr.LinkIndex(i)).L/2;0;0];0 0 0 1]; 
            dx      = Xr(2)-Xr(1);
            Xpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
            Ypatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
            Zpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
            i_patch = 1;

            if Tr.VLinks(Tr.LinkIndex(i)).CS=='C'

                r_fn  = Tr.VLinks(Tr.LinkIndex(i)).r;
                r     = r_fn(0);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = r*sin(theta);
                z     = r*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];

            elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='R'

                h_fn  = Tr.VLinks(Tr.LinkIndex(i)).h;
                w_fn  = Tr.VLinks(Tr.LinkIndex(i)).w;
                h     = h_fn(0);
                w     = w_fn(0);
                x     = [0 0 0 0 0];
                y     = [h/2 -h/2 -h/2 h/2 h/2];
                z     = [w/2 w/2 -w/2 -w/2 w/2];
                pos   = [x;y;z;ones(1,5)];

            elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='E'

                a_fn  = Tr.VLinks(Tr.LinkIndex(i)).a;
                b_fn  = Tr.VLinks(Tr.LinkIndex(i)).b;
                a     = a_fn(0);
                b     = b_fn(0);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = a*sin(theta);
                z     = b*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];
            end

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

                if Tr.VLinks(Tr.LinkIndex(i)).CS=='C'

                    r     = r_fn(Xr(ii)/L);
                    theta = linspace(0,2*pi,n_r);
                    x     = zeros(1,n_r);
                    y     = r*sin(theta);
                    z     = r*cos(theta);
                    pos   = [x;y;z;ones(1,n_r)];

                elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='R'

                    h     = h_fn(Xr(ii)/L);
                    w     = w_fn(Xr(ii)/L);
                    x     = [0 0 0 0 0];
                    y     = [h/2 -h/2 -h/2 h/2 h/2];
                    z     = [w/2 w/2 -w/2 -w/2 w/2];
                    pos   = [x;y;z;ones(1,5)];

                elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='E'

                    a     = h_fn(Xr(ii)/L);
                    b     = w_fn(Xr(ii)/L);
                    theta = linspace(0,2*pi,n_r);
                    x     = zeros(1,n_r);
                    y     = a*sin(theta);
                    z     = b*cos(theta);
                    pos   = [x;y;z;ones(1,n_r)];
                end

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

            patch(Xpatch,Ypatch,Zpatch,color,'EdgeColor','none');
        else
            CustomShapePlot(g_here);
        end
        gf     = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here = g_here*gf;
    end
    end
    dof_start = dof_start+dof_here;
    %=============================================================================
    %Soft link pieces
    for j=1:(Tr.VLinks(Tr.LinkIndex(i)).npie)-1
        
        dof_here   = Tr.CVTwists{i}(j+1).dof;
        Type       = Tr.CVTwists{i}(j+1).Type;
        q_here     = q(dof_start:dof_start+dof_here-1);
        xi_starfn  = Tr.CVTwists{i}(j+1).xi_starfn;
        gi         = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
        Bdof      = Tr.CVTwists{i}(j+1).Bdof;
        Bodr      = Tr.CVTwists{i}(j+1).Bodr;
        Bh     = Tr.CVTwists{i}(j+1).Bh;
        ld     = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        g_here = g_here*gi;

        Xs = linspace(0,1,n_l);
        H  = Xs(2)-Xs(1);
        
        Z = (1/2)*H;          % Zanna quadrature coefficient
        
        Xpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
        Ypatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
        Zpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
        i_patch = 1;
        
        if Tr.VLinks(Tr.LinkIndex(i)).CS == 'C'

            r_fn  = Tr.VLinks(Tr.LinkIndex(i)).r{j};
            r     = r_fn(0);
            theta = linspace(0,2*pi,n_r);
            x     = zeros(1,n_r);
            y     = r*sin(theta);
            z     = r*cos(theta);
            pos   = [x;y;z;ones(1,n_r)];
            
        elseif Tr.VLinks(Tr.LinkIndex(i)).CS == 'R'
            h_fn  = Tr.VLinks(Tr.LinkIndex(i)).h{j};
            w_fn  = Tr.VLinks(Tr.LinkIndex(i)).w{j};
            h     = h_fn(0);
            w     = w_fn(0);
            x     = [0 0 0 0 0];
            y     = [h/2 -h/2 -h/2 h/2 h/2];
            z     = [w/2 w/2 -w/2 -w/2 w/2];
            pos   = [x;y;z;ones(1,5)];
            
        elseif Tr.VLinks(Tr.LinkIndex(i)).CS == 'E'
            a_fn  = Tr.VLinks(Tr.LinkIndex(i)).a{j};
            b_fn  = Tr.VLinks(Tr.LinkIndex(i)).b{j};
            a     = a_fn(0);
            b     = b_fn(0);
            theta = linspace(0,2*pi,n_r);
            x     = zeros(1,n_r);
            y     = a*sin(theta);
            z     = b*cos(theta);
            pos   = [x;y;z;ones(1,n_r)];
        end
        
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

        Lscale = ld;
        
        for ii=1:n_l-1
            if Tr.VLinks(Tr.LinkIndex(i)).CS == 'C'

                r     = r_fn(Xs(ii+1));
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = r*sin(theta);
                z     = r*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];
                
            elseif Tr.VLinks(Tr.LinkIndex(i)).CS == 'R'

                h     = h_fn(Xs(ii+1));
                w     = w_fn(Xs(ii+1));
                x     = [0 0 0 0 0];
                y     = [h/2 -h/2 -h/2 h/2 h/2];
                z     = [w/2 w/2 -w/2 -w/2 w/2];
                pos   = [x;y;z;ones(1,5)];
                
            elseif Tr.VLinks(Tr.LinkIndex(i)).CS == 'E'

                a     = a_fn(Xs(ii+1));
                b     = b_fn(Xs(ii+1));
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = a*sin(theta);
                z     = b*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];
            end
            
            X = Xs(ii);
            
            X_Z = X+Z;

            xi_Zhere  = xi_starfn(X_Z);
            xi_Zhere(1:3) = xi_Zhere(1:3)*ld;
            
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

        patch(Xpatch,Ypatch,Zpatch,color,'EdgeColor','none');
        
        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gf     = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
        g_here = g_here*gf;
        dof_start  = dof_start+dof_here;
    end
    g_Ltip((i-1)*4+1:i*4,:) = g_here;
end
drawnow
end