%Function for the plot of dynamic simulation
%Last modified by Anup Teejo Mathew - 25/05/2021
function plotqqd(S,t,qqd)
close all


PlottingParameters = S.PlotParameters;
% Forward Kinematics: (Product of exponentials)
N    = S.N;

tic
tmax        = max(t);
v           = VideoWriter('.\Dynamics');
FrameRate   = PlottingParameters.FrameRateValue;
v.FrameRate = FrameRate;
open(v);
figure('units','normalized','outerposition',[0 0 1 1]);


for tt=0:1/FrameRate:tmax

    clf
    
    set(gca,'CameraPosition',PlottingParameters.CameraPosition,...
            'CameraTarget',PlottingParameters.CameraTarget,...
            'CameraUpVector',PlottingParameters.CameraUpVector)
    camlight(PlottingParameters.Az_light,PlottingParameters.El_light)

    axis equal
    grid on
    hold on
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')       
    title(strcat('t= ',num2str(tt)))
    axis ([PlottingParameters.X_lim PlottingParameters.Y_lim PlottingParameters.Z_lim]);
    
    
    qqdtt = interp1(t,qqd,tt);
    q     = qqdtt(1:S.ndof)';
    
    f         = 1;
    dof_start = 1;
    g_here    = S.g_ini;
    for i=1:N % number of links
        %joint
        dof_here   = S.Vtwists(f).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        B_here     = S.Vtwists(f).B;
        xi_star    = S.Vtwists(f).xi_star;

        if dof_here==0 %fixed joint (N)
            g_joint    = eye(4);
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
        g_here     = g_here*g_joint;
        
        n_r   = S.VLinks(i).n_r;
        if S.VLinks(i).CS=='R'
            n_r=5;
        end
        n_l   = S.VLinks(i).n_l;
        color = S.VLinks(i).color;
        
        if S.VLinks(i).linktype=='r'
            L          = S.VLinks(i).L;
            g0f        = S.g0{f};
            g_here     = g_here*g0f;
            xR         = linspace(0,L,n_l);
            g_hereR    = g_here*[eye(3) [-S.VLinks(i).cx;0;0];0 0 0 1]; 
            dx         = xR(2)-xR(1);
            Xpatch  = zeros(n_r,n_l);
            Ypatch  = zeros(n_r,n_l);
            Zpatch  = zeros(n_r,n_l);
            i_patch = 1;

            if S.VLinks(i).CS=='C'

                r_fn  = S.VLinks(i).r;
                r     = r_fn(0);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = r*sin(theta);
                z     = r*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];

            elseif S.VLinks(i).CS=='R'

                h_fn  = S.VLinks(i).h;
                w_fn  = S.VLinks(i).w;
                h     = h_fn(0);
                w     = w_fn(0);
                x     = [0 0 0 0 0];
                y     = [h/2 -h/2 -h/2 h/2 h/2];
                z     = [w/2 w/2 -w/2 -w/2 w/2];
                pos   = [x;y;z;ones(1,5)];
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
                
                if S.VLinks(i).CS=='C'

                    r     = r_fn(xR(ii)/L);
                    theta = linspace(0,2*pi,n_r);
                    x     = zeros(1,n_r);
                    y     = r*sin(theta);
                    z     = r*cos(theta);
                    pos   = [x;y;z;ones(1,n_r)];
                    
                elseif S.VLinks(i).CS=='R'

                    h     = h_fn(xR(ii)/L);
                    w     = w_fn(xR(ii)/L);
                    x     = [0 0 0 0 0];
                    y     = [h/2 -h/2 -h/2 h/2 h/2];
                    z     = [w/2 w/2 -w/2 -w/2 w/2];
                    pos   = [x;y;z;ones(1,5)];
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

            g0f(2:3,4) = -g0f(2:3,4);
            g0f(1,4)   = S.VLinks(i).L-g0f(1,4);
            g_here     = g_here*g0f;
            patch(Xpatch,Ypatch,Zpatch,color,'EdgeColor','none');

        end
        
        
        f=f+1;
        dof_start = dof_start+dof_here;
            %=============================================================================
        for j=1:(S.VLinks(i).npie)-1 % for each piece
            
            dof_here   = S.Vtwists(f).dof;
            q_here     = q(dof_start:dof_start+dof_here-1);
            xi_starfn  = S.Vtwists(f).xi_starfn;
            g0f        = S.g0{f};
            Bdof       = S.Vtwists(f).Bdof;
            Bodr       = S.Vtwists(f).Bodr;
            lpf        = S.VLinks(i).lp{j};
            g_here     = g_here*g0f;
               
            xS          = linspace(0,lpf,n_l);
            color       = S.VLinks(i).color;
            H           = xS(2)-xS(1);
            Z1          = (1/2-sqrt(3)/6)*H;          % Zanna quadrature coefficient
            Z2          = (1/2+sqrt(3)/6)*H;          % Zanna quadrature coefficient
            B_Z1here    = zeros(6,dof_here);
            B_Z2here    = zeros(6,dof_here);

            Xpatch  = zeros(n_r,n_l);
            Ypatch  = zeros(n_r,n_l);
            Zpatch  = zeros(n_r,n_l);
            i_patch = 1;
            
            if S.VLinks(i).CS=='C'
                
                r_fn  = S.VLinks(i).r{j};
                r     = r_fn(0);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = r*sin(theta);
                z     = r*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];
            elseif S.VLinks(i).CS=='R'
                h_fn = S.VLinks(i).h{j};
                w_fn = S.VLinks(i).w{j};
                h    = h_fn(0);
                w    = w_fn(0);
                x    = [0 0 0 0 0];
                y    = [h/2 -h/2 -h/2 h/2 h/2];
                z    = [w/2 w/2 -w/2 -w/2 w/2];
                pos  = [x;y;z;ones(1,5)];
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
            
            for ii=1:n_l-1
                
                if S.VLinks(i).CS=='C'
                    r     = r_fn(xS(ii)/lpf);
                    theta = linspace(0,2*pi,n_r);
                    x   = zeros(1,n_r);
                    y   = r*sin(theta);
                    z   = r*cos(theta);
                    pos = [x;y;z;ones(1,n_r)];
                elseif S.VLinks(i).CS=='R'

                    h   = h_fn(xS(ii)/lpf);
                    w   = w_fn(xS(ii)/lpf);
                    x   = [0 0 0 0 0];
                    y   = [h/2 -h/2 -h/2 h/2 h/2];
                    z   = [w/2 w/2 -w/2 -w/2 w/2];
                    pos = [x;y;z;ones(1,5)];
                end
                
                x    = xS(ii);
                x_Z1 = x+Z1;
                x_Z2 = x+Z2;
                
                for jj=1:6
                    for k=1:Bdof(jj)*Bodr(jj)+Bdof(jj)
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
                
                Gamma_here    = (H/2)*(xi_Z1here+xi_Z2here)+...
                                ((sqrt(3)*H^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;
                k_here        = Gamma_here(1:3);
                theta_here    = norm(k_here);
                gh            = variable_expmap(theta_here,Gamma_here);
                g_here        = g_here*gh;

                pos_here = g_here*pos;
                x_here   = pos_here(1,:);
                y_here   = pos_here(2,:);
                z_here   = pos_here(3,:);


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
            g0f(2:3,4)    = -g0f(2:3,4);
            g_here        = g_here*g0f;
            
            dof_start = dof_start+dof_here;
            f=f+1;
        end
    end
    

    drawnow

    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
implay('.\Dynamics.avi')