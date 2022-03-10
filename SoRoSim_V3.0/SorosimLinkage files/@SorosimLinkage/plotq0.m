%Function that plots the reference configuration
%Last modified by Anup Teejo Mathew 02.03.2022

function plotq0(Tr,Lh,Dh,CLh)

if nargin==1
    Lh=1e3;
end

close all

N         = Tr.N;
Lscale    = Tr.PlotParameters.Lscale;% later
VLinks    = Tr.VLinks;
LinkIndex = Tr.LinkIndex;
iLpre     = Tr.iLpre;
ndof      = Tr.ndof;
g_ini     = Tr.g_ini;
CVTwists  = Tr.CVTwists;

view(45,45)
axis equal
grid on
hold on
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

%Forward kinematics
q         = zeros(ndof,1);
dof_start = 1;
g_Ltip    = repmat(eye(4),N,1);

n_sigplot=N;
for i=1:N
    if VLinks(LinkIndex(i)).linktype=='r'
        n_sigplot=n_sigplot+1;
    end
    for j=VLinks(LinkIndex(i)).npie-1
        n_sigplot=n_sigplot+VLinks(LinkIndex(i)).n_l;
    end
end

g         = zeros(4*n_sigplot,4);
i_sigplot     = 1;

scale = 1.1*Lscale;

o0 = [0,0,0]; %spacial frame origin
oX = [scale/5,0,0];oY=[0 scale/5 0];oZ=[0 0 scale/5];

arrow3(o0,oX,'_r-1',scale);
arrow3(o0,oY,'_g-1',scale);
arrow3(o0,oZ,'_b-1',scale);

for i = 1:N

    if iLpre(i)>0
        g_here=g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
    else
        g_here=g_ini((i-1)*4+1:i*4,:);
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

    g((i_sigplot-1)*4+1:i_sigplot*4,:) = g_here;
    i_sigplot                          = i_sigplot+1;

    LineWidthValue = 1;
    color   = VLinks(LinkIndex(i)).color;
%%%%%%%%%%%%%%%%%%%%%%%%%%To Highlight a Link%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if any(Lh==i)
        LineWidthValue = 2;
        if color=='y'
            color='r';
        else
            color='y';
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%Estimate the scale factor for symbols%%%%%%%%%%%%%%%%%%%%
    if VLinks(LinkIndex(i)).CS=='C'
        if VLinks(LinkIndex(i)).linktype=='s'
            Ai = VLinks(LinkIndex(i)).r{1}(0)^2;
        else
            Ai = VLinks(LinkIndex(i)).r(0)^2;
        end
    elseif VLinks(LinkIndex(i)).CS=='R'
        if VLinks(LinkIndex(i)).linktype=='s'
            Ai = VLinks(LinkIndex(i)).w{1}(0)*VLinks(LinkIndex(i)).h{1}(0);
        else
            Ai = VLinks(LinkIndex(i)).w(0)*VLinks(LinkIndex(i)).h(0);
        end
    elseif VLinks(LinkIndex(i)).CS=='E'
        if VLinks(LinkIndex(i)).linktype=='s'
            Ai = VLinks(LinkIndex(i)).a{1}(0)*VLinks(LinkIndex(i)).b{1}(0);
        else
            Ai = VLinks(LinkIndex(i)).a(0)*VLinks(LinkIndex(i)).b(0);
        end
    end
    scale    = (VLinks(LinkIndex(i)).L*Ai)^(1/3); %change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(B_here)

        if any(B_here(1,:))

            PosR     = RotationSymbol*scale;
            PosR     = [PosR ones(size(PosR,1),1)];
            PosR     = g_here*PosR';
            x_here   = PosR(1,:);
            y_here   = PosR(2,:);
            z_here   = PosR(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(B_here(2,:))

            PosR     = RotationSymbol*scale;
            PosR     = [PosR ones(size(PosR,1),1)];
            PosR     = [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*PosR';
            PosR     = g_here*PosR;
            x_here   = PosR(1,:);
            y_here   = PosR(2,:);
            z_here   = PosR(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(B_here(3,:))

            PosR     = RotationSymbol*scale;
            PosR     = [PosR ones(size(PosR,1),1)];
            PosR     = [0 0 1 0;0 1 0 0;-1 0 0 0;0 0 0 1]*PosR';
            PosR     = g_here*PosR;
            x_here   = PosR(1,:);
            y_here   = PosR(2,:);
            z_here   = PosR(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(B_here(4,:))

            PosT     = TranslationSymbol*scale;
            PosT     = [PosT ones(size(PosT,1),1)];
            PosT     = g_here*PosT';
            x_here   = PosT(1,:);
            y_here   = PosT(2,:);
            z_here   = PosT(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(B_here(5,:))

            PosT     = TranslationSymbol*scale;
            PosT     = [PosT ones(size(PosT,1),1)];
            PosT     = [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*PosT';
            PosT     = g_here*PosT;
            x_here   = PosT(1,:);
            y_here   = PosT(2,:);
            z_here   = PosT(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(B_here(6,:))

            PosT     = TranslationSymbol*scale;
            PosT     = [PosT ones(size(PosT,1),1)];
            PosT     = [0 0 1 0;0 1 0 0;-1 0 0 0;0 0 0 1]*PosT';
            PosT     = g_here*PosT;
            x_here   = PosT(1,:);
            y_here   = PosT(2,:);
            z_here   = PosT(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end
    else
        PosF     = FixedSymbol*scale;
        PosF     = [PosF ones(size(PosF,1),1)];
        PosF     = g_here*PosF';
        x_here   = PosF(1,:);
        y_here   = PosF(2,:);
        z_here   = PosF(3,:);
        plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if VLinks(LinkIndex(i)).linktype=='r'
        L       = VLinks(LinkIndex(i)).L;
        gi      = VLinks(LinkIndex(i)).gi;
        g_here  = g_here*gi;

        text(g_here(1,4),g_here(2,4),g_here(3,4),['Link' num2str(i)],'FontSize',10,'FontWeight', 'bold' ,'Color','k','HorizontalAlignment', 'center')

        g((i_sigplot-1)*4+1:i_sigplot*4,:) = g_here;
        i_sigplot                          = i_sigplot+1;

        n_l     = VLinks(LinkIndex(i)).n_l;
        if rem(n_l,2)==0
            n_l=n_l+1;
        end
        xR      = linspace(0,L,n_l);
        g_hereR = g_here*[eye(3) [-VLinks(LinkIndex(i)).gi(1,4);0;0];0 0 0 1];
        dx      = xR(2)-xR(1);

        for ii=1:n_l
            if VLinks(LinkIndex(i)).CS=='C'
                n_r   = VLinks(LinkIndex(i)).n_r;
                r_fn  = VLinks(LinkIndex(i)).r;

                r     = r_fn(xR(ii)/L);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = r*sin(theta);
                z     = r*cos(theta);
                pos_here   = [x;y;z;ones(1,n_r)];

            elseif VLinks(LinkIndex(i)).CS == 'R'
                h_fn  = VLinks(LinkIndex(i)).h;
                w_fn  = VLinks(LinkIndex(i)).w;
                h     = h_fn(xR(ii)/L);
                w     = w_fn(xR(ii)/L);
                x     = [0 0 0 0 0];
                y     = [h/2 -h/2 -h/2 h/2 h/2];
                z     = [w/2 w/2 -w/2 -w/2 w/2];
                pos_here   = [x;y;z;ones(1,5)];
            elseif VLinks(LinkIndex(i)).CS == 'E'
                n_r   = VLinks(LinkIndex(i)).n_r;
                a_fn  = VLinks(LinkIndex(i)).a;
                b_fn  = VLinks(LinkIndex(i)).b;

                a     = a_fn(xR(ii)/L);
                b     = b_fn(xR(ii)/L);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = a*sin(theta);
                z     = b*cos(theta);
                pos_here   = [x;y;z;ones(1,n_r)];
            end

            pos_here = g_hereR*pos_here;
            x_here   = pos_here(1,:);
            y_here   = pos_here(2,:);
            z_here   = pos_here(3,:);

            %plotting rigid link
            plot3(x_here,y_here,z_here,'color',color,'LineWidth',LineWidthValue)


            %Rigid Link local frame
            if ii==(n_l+1)/2

                oX=0.9*[scale,0,0];oY=0.9*[0 scale 0];oZ=0.9*[0 0 scale];

                o0_here=(g_here*[o0';1])';
                oX_here=(g_here*[oX';1])';
                oY_here=(g_here*[oY';1])';
                oZ_here=(g_here*[oZ';1])';

                arrow3(o0_here(1:3),oX_here(1:3),'r-1',scale*5);
                arrow3(o0_here(1:3),oY_here(1:3),'g-1',scale*5);
                arrow3(o0_here(1:3),oZ_here(1:3),'b-1',scale*5);

            end
            g_hereR  = g_hereR*[eye(3) [dx;0;0];0 0 0 1];
        end


        gf     = VLinks(LinkIndex(i)).gf;
        g_here = g_here*gf;
    end

    dof_start = dof_start+dof_here;
%==========================================================================
    %Soft pieces
    if nargin==2
        Dh=1:(VLinks(LinkIndex(i)).npie)-1;
    elseif nargin==4
        Dh = (VLinks(LinkIndex(i)).npie)-1;
    end
    for j = 1:(VLinks(LinkIndex(i)).npie)-1

        if any(Lh==i)&&any(Dh==j)
            LineWidthValue = 2;
            if VLinks(LinkIndex(i)).color=='y'
                color='r';
            else
                color='y';
            end
        else
            LineWidthValue = 1;
            color          = VLinks(LinkIndex(i)).color;
        end

        xi_starfn = VTwists(j+1).xi_starfn;
        gid       = VLinks(LinkIndex(i)).gi{j};
        Bdof      = VTwists(j+1).Bdof;
        Bodr      = VTwists(j+1).Bodr;
        ld        = VLinks(LinkIndex(i)).lp{j};
        g_here    = g_here*gid;
        dof_here  = VTwists(j+1).dof;
        q_here    = q(dof_start:dof_start+dof_here-1);

        g((i_sigplot-1)*4+1:i_sigplot*4,:) = g_here;
        i_sigplot                          = i_sigplot+1;

        %local frame at begining of piece
        oX=0.9*[scale,0,0];oY=0.9*[0 scale 0];oZ=0.9*[0 0 scale];

        o0_here=(g_here*[o0';1])';
        oX_here=(g_here*[oX';1])';
        oY_here=(g_here*[oY';1])';
        oZ_here=(g_here*[oZ';1])';

        arrow3(o0_here(1:3),oX_here(1:3),'r-1',scale*5);
        arrow3(o0_here(1:3),oY_here(1:3),'g-1',scale*5);
        arrow3(o0_here(1:3),oZ_here(1:3),'b-1',scale*5);

        n_l      = VLinks(LinkIndex(i)).n_l;
        Xs       = linspace(0,ld,n_l);

        H        = Xs(2)-Xs(1);
        Z1       = (1/2-sqrt(3)/6)*H;      %Zanna quadrature coefficient
        Z2       = (1/2+sqrt(3)/6)*H;      %Zanna quadrature coefficient
        B_Z1here = zeros(6,dof_here);
        B_Z2here = zeros(6,dof_here);

        if VLinks(LinkIndex(i)).CS == 'C'
            n_r   = VLinks(LinkIndex(i)).n_r;
            r_fn  = VLinks(LinkIndex(i)).r{j};
            r     = r_fn(0);
            theta = linspace(0,2*pi,n_r);
            x     = zeros(1,n_r);
            y     = r*sin(theta);
            z     = r*cos(theta);
            pos_here   = [x;y;z;ones(1,n_r)];
        elseif VLinks(LinkIndex(i)).CS =='R'
            h_fn  = VLinks(LinkIndex(i)).h{j};
            w_fn  = VLinks(LinkIndex(i)).w{j};
            h     = h_fn(0);
            w     = w_fn(0);
            x     = [0 0 0 0 0];
            y     = [h/2 -h/2 -h/2 h/2 h/2];
            z     = [w/2 w/2 -w/2 -w/2 w/2];
            pos_here   = [x;y;z;ones(1,5)];
        elseif VLinks(LinkIndex(i)).CS =='E'
            n_r   = VLinks(LinkIndex(i)).n_r;
            a_fn  = VLinks(LinkIndex(i)).a{j};
            b_fn  = VLinks(LinkIndex(i)).b{j};
            a     = a_fn(0);
            b     = b_fn(0);
            theta = linspace(0,2*pi,n_r);
            x     = zeros(1,n_r);
            y     = a*sin(theta);
            z     = b*cos(theta);
            pos_here   = [x;y;z;ones(1,n_r)];
        end

        pos_here = g_here*pos_here;
        x_here   = pos_here(1,:);
        y_here   = pos_here(2,:);
        z_here   = pos_here(3,:);
        plot3(x_here,y_here,z_here,'color',color,'LineWidth',LineWidthValue)

        divisionnumbernotadded=1;
        for ii = 1:n_l-1

            if VLinks(LinkIndex(i)).CS == 'C'
                r     = r_fn(Xs(ii)/ld);
                %x and theta are defined before
                y     = r*sin(theta);
                z     = r*cos(theta);
                pos_here   = [x;y;z;ones(1,n_r)];
            elseif VLinks(LinkIndex(i)).CS=='R'
                h     = h_fn(Xs(ii)/ld);
                w     = w_fn(Xs(ii)/ld);
                %x is defined before
                y     = [h/2 -h/2 -h/2 h/2 h/2];
                z     = [w/2 w/2 -w/2 -w/2 w/2];
                pos_here   = [x;y;z;ones(1,5)];
            elseif VLinks(LinkIndex(i)).CS=='E'
                a     = a_fn(Xs(ii)/ld);
                b     = b_fn(Xs(ii)/ld);
                %x and theta are defined before
                y     = a*sin(theta);
                z     = b*cos(theta);
                pos_here   = [x;y;z;ones(1,n_r)];
            end

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
                xi_Z1here  = B_Z1here*q_here+xi_starfn(X_Z1/ld);
                xi_Z2here  = B_Z2here*q_here+xi_starfn(X_Z2/ld);
            else
                xi_Z1here  = xi_starfn(X_Z1/ld);
                xi_Z2here  = xi_starfn(X_Z2/ld);
            end

            Gamma_here = (H/2)*(xi_Z1here+xi_Z2here)+...
                ((sqrt(3)*H^2)/12)*dinamico_adj(xi_Z1here)*xi_Z2here;

            k_here     = Gamma_here(1:3);
            theta_here = norm(k_here);
            gh         = variable_expmap(theta_here,Gamma_here);
            g_here     = g_here*gh;

            g((i_sigplot-1)*4+1:i_sigplot*4,:) = g_here;
            i_sigplot                          = i_sigplot+1;

            pos_here   = g_here*pos_here;
            x_here     = pos_here(1,:);
            y_here     = pos_here(2,:);
            z_here     = pos_here(3,:);

            %plotting soft link pieces
            plot3(x_here,y_here,z_here,'color',color,'LineWidth',LineWidthValue)
            %hold on

            if ii>n_l/2&&divisionnumbernotadded
                text(g_here(1,4),g_here(2,4),g_here(3,4),num2str(j),'FontSize',10,'FontWeight', 'bold','Color','k','HorizontalAlignment', 'center')
                divisionnumbernotadded=0;
            end

        end

        %Local frame at the end of piece
        oX=0.9*[scale,0,0];oY=0.9*[0 scale 0];oZ=0.9*[0 0 scale];

        o0_here=(g_here*[o0';1])';
        oX_here=(g_here*[oX';1])';
        oY_here=(g_here*[oY';1])';
        oZ_here=(g_here*[oZ';1])';

        arrow3(o0_here(1:3),oX_here(1:3),'r-1',scale*5);
        arrow3(o0_here(1:3),oY_here(1:3),'g-1',scale*5);
        arrow3(o0_here(1:3),oZ_here(1:3),'b-1',scale*5);

        %updating g, Jacobian, Jacobian_dot and eta at X=L
        gfd    = VLinks(LinkIndex(i)).gf{j};
        g_here = g_here*gfd;

        dof_start  = dof_start+dof_here;
        if j==floor(VLinks(LinkIndex(i)).npie/2)
            text(g_here(1,4),g_here(2,4),g_here(3,4),['Link',num2str(i)],'FontSize',10,'FontWeight','bold','Color','k','HorizontalAlignment', 'center')
        end
    end
    g_Ltip((i-1)*4+1:i*4,:) = g_here;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Adjusting axis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axis tight
xlim([-0.1*max(abs(xlim)),1.1*max(abs(xlim))])
ylim([-max(abs(ylim)),max(abs(ylim))])
zlim([-max(abs(zlim)),max(abs(zlim))])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for iCL=1:Tr.nCLj

    if Tr.iACL(iCL)>0
        g_here=g_Ltip((Tr.iACL(iCL)-1)*4+1:Tr.iACL(iCL)*4,:)*Tr.gACLj{iCL};
    else
        g_here=Tr.gACLj{iCL};
    end

    LineWidthValue = 1;
    if nargin==1||nargin==2
        CLh=1e3;
    end
    if any(CLh==iCL)
        LineWidthValue = 2;
    end
    B_here   = Tr.VTwistsCLj(iCL).B;
    if ~isempty(B_here)

        if any(B_here(1,:))

            PosR     = RotationSymbol*scale;
            PosR     = [PosR ones(size(PosR,1),1)];
            PosR     = g_here*PosR';
            x_here   = PosR(1,:);
            y_here   = PosR(2,:);
            z_here   = PosR(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(B_here(2,:))

            PosR     = RotationSymbol*scale;
            PosR     = [PosR ones(size(PosR,1),1)];
            PosR     = [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*PosR';
            PosR     = g_here*PosR;
            x_here   = PosR(1,:);
            y_here   = PosR(2,:);
            z_here   = PosR(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(B_here(3,:))

            PosR     = RotationSymbol*scale;
            PosR     = [PosR ones(size(PosR,1),1)];
            PosR     = [0 0 1 0;0 1 0 0;-1 0 0 0;0 0 0 1]*PosR';
            PosR     = g_here*PosR;
            x_here   = PosR(1,:);
            y_here   = PosR(2,:);
            z_here   = PosR(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(B_here(4,:))

            PosT     = TranslationSymbol*scale;
            PosT     = [PosT ones(size(PosT,1),1)];
            PosT     = g_here*PosT';
            x_here   = PosT(1,:);
            y_here   = PosT(2,:);
            z_here   = PosT(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(B_here(5,:))

            PosT     = TranslationSymbol*scale;
            PosT     = [PosT ones(size(PosT,1),1)];
            PosT     = [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*PosT';
            PosT     = g_here*PosT;
            x_here   = PosT(1,:);
            y_here   = PosT(2,:);
            z_here   = PosT(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(B_here(6,:))

            PosT     = TranslationSymbol*scale;
            PosT     = [PosT ones(size(PosT,1),1)];
            PosT     = [0 0 1 0;0 1 0 0;-1 0 0 0;0 0 0 1]*PosT';
            PosT     = g_here*PosT;
            x_here   = PosT(1,:);
            y_here   = PosT(2,:);
            z_here   = PosT(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end
    else
        PosF     = FixedSymbol*scale;
        PosF     = [PosF ones(size(PosF,1),1)];
        PosF     = g_here*PosF';
        x_here   = PosF(1,:);
        y_here   = PosF(2,:);
        z_here   = PosF(3,:);
        plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
    end

end


% scale = scale_spacial;
%Gravity direction if present
scale = scale*2;
if Tr.Gravity

    Gdir = (Tr.G(4:6))/norm(Tr.G(4:6));
    g1=[mean(xlim),max(ylim)*0.9,max(zlim)*0.9];
    g2=g1+Gdir'*scale;
    arrow3(g1,g2);

    txt = '   g';
    text(g1(1),g1(2),g1(3),txt,'FontSize',10,'FontWeight','bold')

end

%Point Force direction if present
if Tr.PointForce
    for ii=1:Tr.np

        Fp_vec = Tr.Fp_vec{ii};
        Fp_loc = Tr.Fp_loc{ii};
        Fp_vec = Fp_vec(0); %at time t=0;
        i_sigplot  = 1;

        for i=1:N

            i_sigplot=i_sigplot+1;


            if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
                if i==Fp_loc(1)
                    g_here = g((i_sigplot-1)*4+1:i_sigplot*4,:);

                    if ~Tr.FollowerForce{ii}
                        pos_here = g_here(1:3,4);
                        g_here(1:3,4)=zeros(3,1);
                        g_adj  = dinamico_Adjoint(ginv(g_here));
                        Fp_vec = g_adj*Fp_vec;
                        g_here(1:3,4)=pos_here;
                    end

                    p1     = g_here*[0 0 0 1]';
                    Fdir   = Fp_vec(4:6)/norm(Fp_vec(4:6));
                    Fdir   = g_here*[Fdir;0];
                    Fdir   = Fdir(1:3);
                    Mdir   = Fp_vec(1:3)/norm(Fp_vec(1:3));
                    Mdir   = g_here*[Mdir;0];
                    Mdir   = Mdir(1:3);
                    break;
                end
                i_sigplot = i_sigplot+1;
            end

            for j=1:VLinks(LinkIndex(i)).npie-1

                i_sigplot = i_sigplot+VLinks(LinkIndex(i)).n_l-1;

                if Fp_loc(1)==i&&Fp_loc(2)==j
                    g_here = g((i_sigplot-1)*4+1:i_sigplot*4,:);

                    if ~Tr.FollowerForce{ii}
                        pos_here = g_here(1:3,4);
                        g_here(1:3,4)=zeros(3,1);
                        g_adj  = dinamico_Adjoint(ginv(g_here));
                        Fp_vec = g_adj*Fp_vec;
                        g_here(1:3,4)=pos_here;
                    end

                    p1     = g_here*[0 0 0 1]';
                    Fdir   = Fp_vec(4:6)/norm(Fp_vec(4:6));
                    Fdir   = g_here*[Fdir;0];
                    Fdir   = Fdir(1:3);
                    Mdir   = Fp_vec(1:3)/norm(Fp_vec(1:3));
                    Mdir   = g_here*[Mdir;0];
                    Mdir   = Mdir(1:3);
                end

                i_sigplot = i_sigplot+1;
            end
        end
        p1 = p1(1:3);
        if norm(Fp_vec(4:6))>0
            p2 = p1+Fdir*scale;
            arrow3(p1',p2','f');

            txt = ['   F',num2str(ii)];
            text(p2(1),p2(2),p2(3),txt)
        end

        if norm(Fp_vec(1:3))>0
            p2 = p1+Mdir*scale;
            arrow3(p1',p2','v');
            p2 = p1+Mdir*scale/1.5;
            arrow3(p1',p2','v');

            txt = ['   M',num2str(ii)];
            text(p2(1),p2(2),p2(3),txt)
        end

    end
end


if Tr.Actuated

    for ii=Tr.i_jact

        i_sig=1;
        dof_start=1;
        for i=1:ii-1
            i_sig=i_sig+1;
            if VLinks(LinkIndex(i)).linktype=='r'
                i_sig = i_sig+1;
            end
            dof_start=dof_start+CVTwists{i}(1).dof;
            for j=1:VLinks(LinkIndex(i)).npie-1
                i_sig = i_sig+VLinks(LinkIndex(i)).nGauss{j};
                dof_start=dof_start+Tr.CVTwists{i}(j+1).dof;
            end
        end
        jpos    = g((i_sig-1)*4+1:i_sig*4,:)*[0 0 0 1]';
        i_jactq = Tr.i_jactq;
        WrenchControlledj=Tr.WrenchControlled(find(i_jactq==dof_start));
        if WrenchControlledj
            txt = '  W';
        else
            txt = '  Q';
        end
        text(jpos(1),jpos(2),jpos(3),txt,'FontSize',10,'Color','k')

    end

    if Tr.n_sact>0
        if ~isempty(Tr.CablePoints)
            CablePoints = Tr.CablePoints;
            for ii=1:Tr.n_sact
                for i=1:Tr.N
                    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='s'
                        y_fn=CablePoints.Cy_fn{ii,i};
                        z_fn=CablePoints.Cz_fn{ii,i};
                        CableDisplay(y_fn,z_fn,ii,i,Tr.Sdiv{ii,i},Tr.Ediv{ii,i},Tr,0)
                    end
                end
            end
        end
    end
end

axis tight
drawnow
end