%Function that plots the reference configuration
%Last modified by Anup Teejo Mathew 29.11.2024
%Modified to also render the annular pneumatic chambers of the
%'am_isupport' cross section (via computeChamberYZ) in the soft-piece loop.

function plotq0(Linkage,Lh,Dh,CLh)

if nargin==1
    Lh=1e4; %I guess no one will make a linkage with 10,000 links. Happy to be wrong
end

close all

N         = Linkage.N;
Lscale    = Linkage.PlotParameters.Lscale;% later
VLinks    = Linkage.VLinks;
LinkIndex = Linkage.LinkIndex;
iLpre     = Linkage.iLpre;
g_ini     = Linkage.g_ini;
CVRods    = Linkage.CVRods;

figure
view(45,45)
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

%Forward kinematics

dof_start = 1;
g_Ltip    = repmat(eye(4),N,1);

n_sigplot=N;
for i=1:N
    if VLinks(LinkIndex(i)).linktype=='r'&&VLinks(LinkIndex(i)).L>0
        n_sigplot=n_sigplot+1;
    end
    for j=1:VLinks(LinkIndex(i)).npie-1
        n_sigplot=n_sigplot+VLinks(LinkIndex(i)).n_l;
    end
end

g         = zeros(4*n_sigplot,4);
i_sigplot = 1;

Lscale(Lscale==0) = 1;
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

    VRods = CVRods{i};

    %Joints
    xi_star  = VRods(1).xi_star;
    Phi_here = VRods(1).Phi; %Joint Phi
    dof_here = VRods(1).dof;

    if dof_here == 0 %fixed joint (N)
        g_joint = eye(4);
    else
        xi         = xi_star;
        g_joint    = variable_expmap_g(xi);
    end

    g_here = g_here*g_joint;

    g((i_sigplot-1)*4+1:i_sigplot*4,:) = g_here;
    i_sigplot                          = i_sigplot+1;

    scale = VLinks(LinkIndex(i)).Lscale;
    
    %plot a frame at X=1 of the joint
    oX=0.9*[scale,0,0];oY=0.9*[0 scale 0];oZ=0.9*[0 0 scale];
    
    o0_here=(g_here*[o0';1])';
    oX_here=(g_here*[oX';1])';
    oY_here=(g_here*[oY';1])';
    oZ_here=(g_here*[oZ';1])';
    
    arrow3(o0_here(1:3),oX_here(1:3),'r-1',scale*5);
    arrow3(o0_here(1:3),oY_here(1:3),'g-1',scale*5);
    arrow3(o0_here(1:3),oZ_here(1:3),'b-1',scale*5);

    LineWidthValue = 1;
    color = VLinks(LinkIndex(i)).color;

    %To Highlight a Link
    if any(Lh==i)
        LineWidthValue = 2;% thicker
        color = [0.7, 1, 0];% lime green
    end

    %% Plotting the Joint
    if ~isempty(Phi_here)

        if any(Phi_here(1,:))

            PosR     = RotationSymbol*scale;
            PosRH    = [PosR ones(size(PosR,1),1)];
            PosRH    = g_here*PosRH';
            x_here   = PosRH(1,:);
            y_here   = PosRH(2,:);
            z_here   = PosRH(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(2,:))

            PosR     = RotationSymbol*scale;
            PosRH    = [PosR ones(size(PosR,1),1)];
            PosRH    = [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*PosRH';
            PosRH    = g_here*PosRH;
            x_here   = PosRH(1,:);
            y_here   = PosRH(2,:);
            z_here   = PosRH(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(3,:))

            PosR     = RotationSymbol*scale;
            PosRH    = [PosR ones(size(PosR,1),1)];
            PosRH    = [0 0 1 0;0 1 0 0;-1 0 0 0;0 0 0 1]*PosRH';
            PosRH    = g_here*PosRH;
            x_here   = PosRH(1,:);
            y_here   = PosRH(2,:);
            z_here   = PosRH(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(4,:))

            PosT     = TranslationSymbol*scale;
            PosTH    = [PosT ones(size(PosT,1),1)];
            PosTH    = g_here*PosTH';
            x_here   = PosTH(1,:);
            y_here   = PosTH(2,:);
            z_here   = PosTH(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(5,:))

            PosT     = TranslationSymbol*scale;
            PosTH    = [PosT ones(size(PosT,1),1)];
            PosTH    = [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*PosTH';
            PosTH     = g_here*PosTH;
            x_here   = PosTH(1,:);
            y_here   = PosTH(2,:);
            z_here   = PosTH(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(6,:))

            PosT     = TranslationSymbol*scale;
            PosTH     = [PosT ones(size(PosT,1),1)];
            PosTH     = [0 0 1 0;0 1 0 0;-1 0 0 0;0 0 0 1]*PosTH';
            PosTH     = g_here*PosTH;
            x_here   = PosTH(1,:);
            y_here   = PosTH(2,:);
            z_here   = PosTH(3,:);
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
        end
    else
        PosF     = FixedSymbol*scale;
        PosFH    = [PosF ones(size(PosF,1),1)];
        PosFH    = g_here*PosFH';
        x_here   = PosFH(1,:);
        y_here   = PosFH(2,:);
        z_here   = PosFH(3,:);
        plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)
    end

%%
    if VLinks(LinkIndex(i)).L>0&&VLinks(LinkIndex(i)).linktype=='r'

        L      = VLinks(LinkIndex(i)).L;
        gi     = VLinks(LinkIndex(i)).gi;
        g_here = g_here*gi;

        text(g_here(1,4),g_here(2,4),g_here(3,4),['Link' num2str(i)],'FontSize',10,'FontWeight', 'bold' ,'Color','k','HorizontalAlignment', 'center')

        g((i_sigplot-1)*4+1:i_sigplot*4,:) = g_here;
        i_sigplot                          = i_sigplot+1;

        n_l     = VLinks(LinkIndex(i)).n_l;
        n_r   = VLinks(LinkIndex(i)).n_r;
        if rem(n_l,2)==0
            n_l=n_l+1;
        end

        if ~Linkage.VLinks(Linkage.LinkIndex(i)).CPF

            Xr      = linspace(0,L,n_l);
            g_hereR = g_here*[eye(3) [-VLinks(LinkIndex(i)).cx;0;0];0 0 0 1];
            dx      = Xr(2)-Xr(1);

            for ii=1:n_l
                [y,z] = computeBoundaryYZ(VLinks(LinkIndex(i)),Xr(ii)/L);
                pos  = [zeros(1,n_r);y;z;ones(1,n_r)]; %homogeneous positions in local frame 4xn_r

                pos = g_hereR*pos;
                x_here   = pos(1,:);
                y_here   = pos(2,:);
                z_here   = pos(3,:);

                %plotting rigid link
                plot3(x_here,y_here,z_here,'color',color,'LineWidth',LineWidthValue)
                
                g_hereR  = g_hereR*[eye(3) [dx;0;0];0 0 0 1];
            end
        else
            CustomShapePlot(g_here);
        end
        %Rigid Link local frame

        oX=0.9*[scale,0,0];oY=0.9*[0 scale 0];oZ=0.9*[0 0 scale];

        o0_here=(g_here*[o0';1])';
        oX_here=(g_here*[oX';1])';
        oY_here=(g_here*[oY';1])';
        oZ_here=(g_here*[oZ';1])';

        arrow3(o0_here(1:3),oX_here(1:3),'r-1',scale*5);
        arrow3(o0_here(1:3),oY_here(1:3),'g-1',scale*5);
        arrow3(o0_here(1:3),oZ_here(1:3),'b-1',scale*5);
        
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

        if any(Lh==i)&&any(Dh==j) %Highlight color
            LineWidthValue = 2;
            color = [0.7, 1, 0];% lime green
        else
            LineWidthValue = 1;
            color          = VLinks(LinkIndex(i)).color;
        end

        xi_starfn = VRods(j+1).xi_starfn;
        gid       = VLinks(LinkIndex(i)).gi{j};

        dof_here  = VRods(j+1).dof;
        g_here    = g_here*gid;

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

        ld      = Linkage.VLinks(Linkage.LinkIndex(i)).ld{j};

        n_l     = VLinks(LinkIndex(i)).n_l;
        n_r     = VLinks(LinkIndex(i)).n_r;

        Xs      = linspace(0,1,n_l);

        H       = Xs(2)-Xs(1);
        Z       = 0.5*H;      %2nd order Zanna quadrature coefficient

        is_amis = strcmp(VLinks(LinkIndex(i)).CS,'am_isupport'); %AM I-Support chambers

        %cross sectional shape Circular, Rectangular, and Ellipsoidal
        [y,z] = computeBoundaryYZ(VLinks(LinkIndex(i)),0,j);
        pos  = [zeros(1,n_r);y;z;ones(1,n_r)]; %homogeneous positions in local frame 4xn_r

        pos = g_here*pos;
        x_here   = pos(1,:);
        y_here   = pos(2,:);
        z_here   = pos(3,:);
        if ~is_amis %skip envelope outline for AM I-Support (chambers drawn instead)
            plot3(x_here,y_here,z_here,'color',color,'LineWidth',LineWidthValue)
        end

        %AM I-Support: draw the annular pneumatic chambers at X=0
        if is_amis
            [Yc,Zc] = computeChamberYZ(VLinks(LinkIndex(i)),0,j);
            for k=1:numel(Yc)
                posc = g_here*[zeros(1,numel(Yc{k}));Yc{k};Zc{k};ones(1,numel(Yc{k}))];
                plot3(posc(1,:),posc(2,:),posc(3,:),'color',color,'LineWidth',LineWidthValue)
            end
        end

        divisionnumbernotadded=true;
        for ii = 1:n_l-1
            %cross sectional shape Circular, Rectangular, and Ellipsoidal
            [y,z] = computeBoundaryYZ(VLinks(LinkIndex(i)),Xs(ii+1),j);
            pos  = [zeros(1,n_r);y;z;ones(1,n_r)]; %homogeneous positions in local frame 4xn_r

            X   = Xs(ii);
            X_Z = X+Z;

            xi_Zhere  = xi_starfn(X_Z);

            Gamma_here = H*ld*xi_Zhere;
            gh         = variable_expmap_g(Gamma_here);
            g_here     = g_here*gh;

            g((i_sigplot-1)*4+1:i_sigplot*4,:) = g_here;
            i_sigplot                          = i_sigplot+1;

            pos   = g_here*pos;
            x_here     = pos(1,:);
            y_here     = pos(2,:);
            z_here     = pos(3,:);

            %plotting soft link pieces
            if ~is_amis %skip envelope outline for AM I-Support (chambers drawn instead)
                plot3(x_here,y_here,z_here,'color',color,'LineWidth',LineWidthValue)
            end
            %hold on

            %AM I-Support: draw the annular pneumatic chambers at this station
            if is_amis
                [Yc,Zc] = computeChamberYZ(VLinks(LinkIndex(i)),Xs(ii+1),j);
                for k=1:numel(Yc)
                    posc = g_here*[zeros(1,numel(Yc{k}));Yc{k};Zc{k};ones(1,numel(Yc{k}))];
                    plot3(posc(1,:),posc(2,:),posc(3,:),'color',color,'LineWidth',LineWidthValue)
                end
            end

            if ii>n_l/2&&divisionnumbernotadded
                text(g_here(1,4),g_here(2,4),g_here(3,4),num2str(j),'FontSize',10,'FontWeight', 'bold','Color','k','HorizontalAlignment', 'center')
                divisionnumbernotadded=false;
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

%% Closed Chain

for iCL=1:Linkage.nCLj
% Frame A
    if Linkage.iCLA(iCL)>0
        g_here=g_Ltip((Linkage.iCLA(iCL)-1)*4+1:Linkage.iCLA(iCL)*4,:)*Linkage.gACLj{iCL};
    else
        g_here=Linkage.gACLj{iCL};
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

    LineWidthValue = 1;
    if nargin==1||nargin==2
        CLh=1e4; % guessing no more than 10,000 closed loop. Happy to be wrong
    end
    if any(CLh==iCL)
        LineWidthValue = 2;
    end

% Plot Closed Loop joint at A

    Phi_here   = Linkage.VRodsCLj(iCL).Phi;
    if ~isempty(Phi_here)

        if any(Phi_here(1,:))

            PosR     = RotationSymbol*scale;
            PosRH     = [PosR ones(size(PosR,1),1)];
            PosRH     = g_here*PosRH';
            x_here   = PosRH(1,:);
            y_here   = PosRH(2,:);
            z_here   = PosRH(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(2,:))

            PosR     = RotationSymbol*scale;
            PosRH     = [PosR ones(size(PosR,1),1)];
            PosRH     = [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*PosRH';
            PosRH     = g_here*PosRH;
            x_here   = PosRH(1,:);
            y_here   = PosRH(2,:);
            z_here   = PosRH(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(3,:))

            PosR     = RotationSymbol*scale;
            PosRH     = [PosR ones(size(PosR,1),1)];
            PosRH     = [0 0 1 0;0 1 0 0;-1 0 0 0;0 0 0 1]*PosRH';
            PosRH     = g_here*PosRH;
            x_here   = PosRH(1,:);
            y_here   = PosRH(2,:);
            z_here   = PosRH(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(4,:))

            PosT     = TranslationSymbol*scale;
            PosTH     = [PosT ones(size(PosT,1),1)];
            PosTH     = g_here*PosTH';
            x_here   = PosTH(1,:);
            y_here   = PosTH(2,:);
            z_here   = PosTH(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(5,:))

            PosT     = TranslationSymbol*scale;
            PosTH     = [PosT ones(size(PosT,1),1)];
            PosTH     = [0 -1 0 0;1 0 0 0;0 0 1 0;0 0 0 1]*PosTH';
            PosTH     = g_here*PosTH;
            x_here   = PosTH(1,:);
            y_here   = PosTH(2,:);
            z_here   = PosTH(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end

        if any(Phi_here(6,:))

            PosT     = TranslationSymbol*scale;
            PosTH     = [PosT ones(size(PosT,1),1)];
            PosTH     = [0 0 1 0;0 1 0 0;-1 0 0 0;0 0 0 1]*PosTH';
            PosTH     = g_here*PosTH;
            x_here   = PosTH(1,:);
            y_here   = PosTH(2,:);
            z_here   = PosTH(3,:);
            plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
        end
    else
        PosF     = FixedSymbol*scale;
        PosFH     = [PosF ones(size(PosF,1),1)];
        PosFH     = g_here*PosFH';
        x_here   = PosFH(1,:);
        y_here   = PosFH(2,:);
        z_here   = PosFH(3,:);
        plot3(x_here,y_here,z_here,'--k','LineWidth',LineWidthValue)
    end

    % Frame B

    if Linkage.iCLB(iCL)>0
        g_here=g_Ltip((Linkage.iCLB(iCL)-1)*4+1:Linkage.iCLB(iCL)*4,:)*Linkage.gBCLj{iCL};
    else
        g_here=Linkage.gBCLj{iCL};
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

end


%%

% scale = scale_spacial;
%Gravity direction if present
scale = scale*2; %dont know why I did this

if Linkage.Gravity

    Gdir = (Linkage.G(4:6))/norm(Linkage.G(4:6));
    axis tight
    g1=[mean(xlim),max(ylim)*0.9,max(zlim)*0.9];
    g2=g1+Gdir'*scale;
    arrow3(g1,g2);

    txt = '   G';
    text(g1(1),g1(2),g1(3),txt,'FontSize',10,'FontWeight','bold')

end

%Point Force direction if present

for ii=1:Linkage.np

    Fp_vec = Linkage.Fp_vec{ii};
    Fp_loc = Linkage.Fp_loc{ii};
    Fp_vec = Fp_vec(0); %at time t=0;
    i_sigplot  = 1;

    for i=1:N

        i_sigplot=i_sigplot+1;

        if Linkage.VLinks(Linkage.LinkIndex(i)).L>0
            if Linkage.VLinks(Linkage.LinkIndex(i)).linktype=='r'
                if i==Fp_loc(1)
                    g_here = g((i_sigplot-1)*4+1:i_sigplot*4,:);

                    if ~Linkage.LocalWrench(ii)
                        pos = g_here(1:3,4);
                        g_here(1:3,4)=zeros(3,1); % only rotating
                        Fp_vec = dinamico_coAdjoint(ginv(g_here))*Fp_vec;
                        g_here(1:3,4)=pos;
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
        end

        for j=1:VLinks(LinkIndex(i)).npie-1

            if Fp_loc(1)==i&&Fp_loc(2)==j
                X = linspace(0,1,n_l);
                for ip=1:VLinks(LinkIndex(i)).n_l
                    if X(ip)>=Fp_loc(3)% if not equal it will be approximate
                        if X(ip)==Fp_loc(3)
                            g_here = g((i_sigplot-1)*4+1:i_sigplot*4,:);
                        else
                            i2 = i_sigplot; i1 = i_sigplot-1;
                            X2 = X(ip); X1 = X(ip-1);
                            Xp = Fp_loc(3);
                            alpha = (X2-Xp)/(X2-X1);
                            g_here = g((i1-1)*4+1:i1*4,:)*alpha+g((i2-1)*4+1:i2*4,:)*(1-alpha);
                        end

                        if ~Linkage.LocalWrench(ii)
                            pos = g_here(1:3,4);
                            g_here(1:3,4)=zeros(3,1);
                            Fp_vec = dinamico_coAdjoint(ginv(g_here))*Fp_vec;
                            g_here(1:3,4)=pos;
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
            else
                i_sigplot = i_sigplot+VLinks(LinkIndex(i)).n_l;
            end
        end
    end

    p1 = p1(1:3);
    if norm(Fp_vec(4:6))>0
        p2 = p1+Fdir*scale;
        arrow3(p1',p2','f');

        txt = ['   f',num2str(ii)];
        text(p2(1),p2(2),p2(3),txt)
    end

    if norm(Fp_vec(1:3))>0
        p2 = p1+Mdir*scale;
        arrow3(p1',p2','v');
        p2 = p1+Mdir*scale/1.5;
        arrow3(p1',p2','v');

        txt = ['   m',num2str(ii)];
        text(p2(1),p2(2),p2(3),txt)
    end

end

%%
if Linkage.Actuated

    for ii=Linkage.i_jact

        i_sig=1;
        dof_start=1;
        for i=1:ii-1
            i_sig=i_sig+1;
            if VLinks(LinkIndex(i)).linktype=='r'
                i_sig = i_sig+1;
            end
            dof_start=dof_start+CVRods{i}(1).dof;
            for j=1:VLinks(LinkIndex(i)).npie-1
                i_sig = i_sig+Linkage.CVRods{i}(j+1).nip;
                dof_start=dof_start+Linkage.CVRods{i}(j+1).dof;
            end
        end
        jpos    = g((i_sig-1)*4+1:i_sig*4,:)*[0 0 0 1]';
        i_jactq = Linkage.i_jactq;
        WrenchControlledj=Linkage.WrenchControlled(find(i_jactq==dof_start));
        if WrenchControlledj
            txt = '  W';
        else
            txt = '  Q';
        end
        text(jpos(1),jpos(2),jpos(3),txt,'FontSize',10,'Color','k')

    end

    if Linkage.n_sact>0
        if ~isempty(Linkage.CableFunction)
            dc_fn = Linkage.CableFunction.dc_fn;
            for ii=1:Linkage.n_sact
                for i=1:Linkage.N
                    if Linkage.VLinks(Linkage.LinkIndex(i)).linktype=='s'
                        dc_fn_here=dc_fn{ii,i};
                        CableDisplay(Linkage,dc_fn_here,i,Linkage.Sdiv(ii,i),Linkage.Ediv(ii,i))
                    end
                end
            end
        end
    end
end

hText = findobj(gca, 'Type', 'Text');
set(hText, 'Interpreter', 'latex');
axis tight
drawnow
end