%Function that plots the reference configuration
%(22.05.2021)

function plotq0(S,ff)
if nargin==1
    ff=1e3;
end
close all
Ltot = S.Ltot;% later

%Plot options
view(45,45)
axis equal
grid on
hold on
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

axis ([-1.1*Ltot 1.1*Ltot...
    -1.1*Ltot 1.1*Ltot...
    -1.1*Ltot 1.1*Ltot])

%Forward kinematics
q         = zeros(S.ndof,1);
N         = S.N;
f         = 1;
dof_start = 1;
g_here    = S.g_ini;

n_sigplot=0;
for i=1:N
    n_sigplot=n_sigplot+1;
    for j=S.VLinks(i).npie-1
        n_sigplot=S.VLinks(i).n_l;
    end
end
g         = zeros(4*n_sigplot,4);
i_sig     = 1;
% ar_l=Ltot/4; %Ha=((ar_l)/un)/3; Wa=Ha/3;

if S.VLinks(1).CS=='C'
    if S.VLinks(1).linktype=='s'
        Ai = S.VLinks(1).r{1}(0)^2;
    else
        Ai = S.VLinks(1).r(0)^2;
    end
else
    if S.VLinks(1).linktype=='s'
        Ai = S.VLinks(1).w{1}(0)*S.VLinks(1).h{1}(0);
    else
        Ai = S.VLinks(1).w(0)*S.VLinks(1).h(0);
    end
end

scale_spacial = 1.1*(S.Ltot*Ai)^(1/3);

scale = scale_spacial;
%Spacial frame origin
o0 = [0,0,0];
oX = [scale,0,0];oY=[0 scale 0];oZ=[0 0 scale];

PBdiag = sqrt(3)*2.2*Ltot;
arrow3(o0,oX,'_r-1',scale/(PBdiag/72)/5);
arrow3(o0,oY,'_g-1',scale/(PBdiag/72)/5);
arrow3(o0,oZ,'_b-1',scale/(PBdiag/72)/5);



% ar_l=Ltot/8;

for i = 1:N
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
    
    g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
    i_sig                         = i_sig+1;
    
    if S.VLinks(i).linktype=='s'
        text(g_here(1,4),g_here(2,4),g_here(3,4),num2str(f),'FontSize',20,'Color','k','HorizontalAlignment', 'center')
    end
    LineWidthValue = 1;
    color   = S.VLinks(i).color;
    if f==ff
        LineWidthValue = 2;
        if color=='y'
            color='r';
        else
            color='y';
        end
    end
    
    if S.VLinks(i).CS=='C'
        if S.VLinks(i).linktype=='s'
            Ai = S.VLinks(i).r{1}(0)^2;
        else
            Ai = S.VLinks(i).r(0)^2;
        end
    else
        if S.VLinks(i).linktype=='s'
            Ai = S.VLinks(i).w{1}(0)*S.VLinks(i).h{1}(0);
        else
            Ai = S.VLinks(i).w(0)*S.VLinks(i).h(0);
        end
    end
    scale    = (S.VLinks(i).L*Ai)^(1/3);
    
    if ~isempty(B_here)
    
        if any(B_here(1,:))

            PosR     = RotationSymbol*scale;
            PosR     = [PosR ones(size(PosR,1),1)];
            PosR     = g_here*PosR';
            x_here   = PosR(1,:);
            y_here   = PosR(2,:);
            z_here   = PosR(3,:);

            %plotting rigid link
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

            %plotting rigid link
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

            %plotting rigid link
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)

        end

        if any(B_here(4,:))

            PosT     = TranslationSymbol*scale;
            PosT     = [PosT ones(size(PosT,1),1)];
            PosT     = g_here*PosT';
            x_here   = PosT(1,:);
            y_here   = PosT(2,:);
            z_here   = PosT(3,:);

            %plotting rigid link
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

            %plotting rigid link
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

            %plotting rigid link
            plot3(x_here,y_here,z_here,'k','LineWidth',LineWidthValue)

        end
    end
    
    if S.VLinks(i).linktype=='r'
        L       = S.VLinks(i).L;
        g0f     = S.g0{f};
        g_here  = g_here*g0f;
        
        g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
        i_sig                         = i_sig+1;
        
        n_l     = S.VLinks(i).n_l;
        if rem(n_l,2)==0
            n_l=n_l+1;
        end
        xR      = linspace(0,L,n_l);
        g_hereR = g_here*[eye(3) [-S.VLinks(i).cx;0;0];0 0 0 1];
        dx      = xR(2)-xR(1);
        
        for ii=1:n_l
            if S.VLinks(i).CS=='C'
                n_r   = S.VLinks(i).n_r;
                r_fn  = S.VLinks(i).r;
                
                r     = r_fn(xR(ii)/L);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = r*sin(theta);
                z     = r*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];
                
            elseif S.VLinks(i).CS == 'R'
                h_fn  = S.VLinks(i).h;
                w_fn  = S.VLinks(i).w;
                h     = h_fn(xR(ii)/L);
                w     = w_fn(xR(ii)/L);
                x     = [0 0 0 0 0];
                y     = [h/2 -h/2 -h/2 h/2 h/2];
                z     = [w/2 w/2 -w/2 -w/2 w/2];
                pos   = [x;y;z;ones(1,5)];
            end
            
            pos_here = g_hereR*pos;
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

                arrow3(o0_here(1:3),oX_here(1:3),'r-1',scale/(PBdiag/72)/5);
                arrow3(o0_here(1:3),oY_here(1:3),'g-1',scale/(PBdiag/72)/5);
                arrow3(o0_here(1:3),oZ_here(1:3),'b-1',scale/(PBdiag/72)/5);
                
                text(g_here(1,4),g_here(2,4),g_here(3,4),num2str(f),'FontSize',20,'Color','k','HorizontalAlignment', 'center')
            end
            g_hereR  = g_hereR*[eye(3) [dx;0;0];0 0 0 1];
        end
        
        g0f(2:3,4) = -g0f(2:3,4);
        g0f(1,4)   = S.VLinks(i).L-g0f(1,4);
        g_here     = g_here*g0f;
    end
    f=f+1;
    dof_start = dof_start+dof_here;
    %=============================================================================
    %Soft pieces
    for j = 1:(S.VLinks(i).npie)-1
        
        if any(ff==f)
            LineWidthValue = 2;
            if S.VLinks(i).color=='y'
                color='r';
            else
                color='y';
            end
        else
            LineWidthValue = 1;
            color          = S.VLinks(i).color;
        end
        
        xi_starfn = S.Vtwists(f).xi_starfn;
        g0f       = S.g0{f};
        Bdof      = S.Vtwists(f).Bdof;
        Bodr      = S.Vtwists(f).Bodr;
        lpf       = S.VLinks(i).lp{j};
        g_here    = g_here*g0f;
        dof_here  = S.Vtwists(f).dof;
        q_here    = q(dof_start:dof_start+dof_here-1);
        
        g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
        i_sig                         = i_sig+1;
        
        n_l      = S.VLinks(i).n_l;
        xS       = linspace(0,lpf,n_l);
        
        H        = xS(2)-xS(1);
        Z1       = (1/2-sqrt(3)/6)*H;      %Zanna quadrature coefficient
        Z2       = (1/2+sqrt(3)/6)*H;      %Zanna quadrature coefficient
        B_Z1here = zeros(6,dof_here);
        B_Z2here = zeros(6,dof_here);
        
        if S.VLinks(i).CS == 'C'
            n_r   = S.VLinks(i).n_r;
            r_fn  = S.VLinks(i).r{j};
            r     = r_fn(0);
            theta = linspace(0,2*pi,n_r);
            x     = zeros(1,n_r);
            y     = r*sin(theta);
            z     = r*cos(theta);
            pos   = [x;y;z;ones(1,n_r)];
            
        elseif S.VLinks(i).CS =='R'
            h_fn  = S.VLinks(i).h{j};
            w_fn  = S.VLinks(i).w{j};
            h     = h_fn(0);
            w     = w_fn(0);
            x     = [0 0 0 0 0];
            y     = [h/2 -h/2 -h/2 h/2 h/2];
            z     = [w/2 w/2 -w/2 -w/2 w/2];
            pos   = [x;y;z;ones(1,5)];
        end
        
        pos_here = g_here*pos;
        x_here   = pos_here(1,:);
        y_here   = pos_here(2,:);
        z_here   = pos_here(3,:);
        plot3(x_here,y_here,z_here,'color',color,'LineWidth',LineWidthValue)
        
        %local frame at begining of piece
        oX=0.9*[scale,0,0];oY=0.9*[0 scale 0];oZ=0.9*[0 0 scale];

        o0_here=(g_here*[o0';1])';
        oX_here=(g_here*[oX';1])';
        oY_here=(g_here*[oY';1])';
        oZ_here=(g_here*[oZ';1])';

        arrow3(o0_here(1:3),oX_here(1:3),'r-1',scale/(PBdiag/72)/5);
        arrow3(o0_here(1:3),oY_here(1:3),'g-1',scale/(PBdiag/72)/5);
        arrow3(o0_here(1:3),oZ_here(1:3),'b-1',scale/(PBdiag/72)/5);
        
        piecenumbernotadded=1;
        for ii = 1:n_l-1
            if S.VLinks(i).CS == 'C'
                n_r   = S.VLinks(i).n_r;
                r     = r_fn(xS(ii)/lpf);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = r*sin(theta);
                z     = r*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];
            elseif S.VLinks(i).CS=='R'
                h     = h_fn(xS(ii)/lpf);
                w     = w_fn(xS(ii)/lpf);
                x     = [0 0 0 0 0];
                y     = [h/2 -h/2 -h/2 h/2 h/2];
                z     = [w/2 w/2 -w/2 -w/2 w/2];
                pos   = [x;y;z;ones(1,5)];
            end
            
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
            
            g((i_sig-1)*4+1:i_sig*4,:)    = g_here;
            i_sig                         = i_sig+1;

            pos_here   = g_here*pos;
            x_here     = pos_here(1,:);
            y_here     = pos_here(2,:);
            z_here     = pos_here(3,:);
            
            %plotting soft link pieces
            plot3(x_here,y_here,z_here,'color',color,'LineWidth',LineWidthValue)
            %hold on
            
            if ii>n_l/2&&piecenumbernotadded
                text(g_here(1,4),g_here(2,4),g_here(3,4),num2str(f),'FontSize',20,'Color','k','HorizontalAlignment', 'center')
                piecenumbernotadded=0;
            end
        end
        
        %Local frame at the end of piece
        oX=0.9*[scale,0,0];oY=0.9*[0 scale 0];oZ=0.9*[0 0 scale];

        o0_here=(g_here*[o0';1])';
        oX_here=(g_here*[oX';1])';
        oY_here=(g_here*[oY';1])';
        oZ_here=(g_here*[oZ';1])';

        arrow3(o0_here(1:3),oX_here(1:3),'r-1',scale/(PBdiag/72)/5);
        arrow3(o0_here(1:3),oY_here(1:3),'g-1',scale/(PBdiag/72)/5);
        arrow3(o0_here(1:3),oZ_here(1:3),'b-1',scale/(PBdiag/72)/5);
        
        %updating g, Jacobian, Jacobian_dot and eta at X=L
        g0f(2:3,4) = -g0f(2:3,4);
        g_here     = g_here*g0f;
        
        dof_start  = dof_start+dof_here;
        LineWidthValue = 1;
        f=f+1;
    end
end

axis tight
xlim([-0.1*max(abs(xlim)),1.1*max(abs(xlim))])
ylim([-max(abs(ylim)),max(abs(ylim))])
zlim([-max(abs(zlim)),max(abs(zlim))])


scale = scale_spacial;
%Gravity direction if present
if S.Gravity
    
    Gdir = (S.G(4:6))/norm(S.G(4:6));
    g1=[mean(xlim),max(ylim)*0.9,max(zlim)*0.9];
    g2=g1+Gdir'*scale;
    arrow3(g1,g2);
    
    txt = '   g';
    text(g1(1),g1(2),g1(3),txt)
    
end

%Gravity direction if present
if S.PointForce
    
    for ii=1:S.np
        
        Fp_vec = S.Fp_vec{ii};
        Fp_loc = S.Fp_loc{ii};
        
        Fp_vec = Fp_vec(0); %at time t=0;
        
        i_sig  = 1;
        f      = 1;

        for i=1:N

            i_sig=i_sig+1;

            if f==Fp_loc
                g_here = g((i_sig-1)*4+1:i_sig*4,:);
                p1     = g_here*[0 0 0 1]';
                Fdir   = Fp_vec(4:6)/norm(Fp_vec(4:6));
                Fdir   = g_here*[Fdir;0];
                Fdir   = Fdir(1:3);
                Mdir   = Fp_vec(1:3)/norm(Fp_vec(1:3));
                Mdir   = g_here*[Mdir;0];
                Mdir   = Mdir(1:3);
            end

            if S.VLinks(i).linktype=='r'
                i_sig = i_sig+1;
            end
            f=f+1;

            if f>Fp_loc
                break;
            end

            for j=1:S.VLinks(i).npie-1

                i_sig = i_sig+S.VLinks(i).n_l-1;

                if f==Fp_loc
                    g_here = g((i_sig-1)*4+1:i_sig*4,:);
                    p1     = g_here*[0 0 0 1]';
                    Fdir   = Fp_vec(4:6)/norm(Fp_vec(4:6));
                    Fdir   = g_here*[Fdir;0];
                    Fdir   = Fdir(1:3);
                    Mdir   = Fp_vec(1:3)/norm(Fp_vec(1:3));
                    Mdir   = g_here*[Mdir;0];
                    Mdir   = Mdir(1:3);
                end

                i_sig = i_sig+1;
                f     = f+1;

                if f>Fp_loc
                    break;
                end

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

            txt = ['   M',num2str(ii)];
            text(p2(1),p2(2),p2(3),txt)
        end

        
    end
    
end

if S.Actuated
    
    for ii=S.i_jact
        
        i_sig=1;
        f=1;
        dof_start=1;
        for i=1:ii-1
            i_sig=i_sig+1;
            if S.VLinks(i).linktype=='r'
                i_sig = i_sig+1;
            end
            dof_start=dof_start+S.Vtwists(f).dof;
            f=f+1;
            for j=1:S.VLinks(i).npie-1
                i_sig = i_sig+S.VLinks(i).nGauss{j};
                dof_start=dof_start+S.Vtwists(f).dof;
                f=f+1;
            end
        end
        jpos    = g((i_sig-1)*4+1:i_sig*4,:)*[0 0 0 1]';
        i_jactq = S.i_jactq;
        WrenchControlledj=S.WrenchControlled(find(i_jactq==dof_start));
        if WrenchControlledj
            txt = '  W';
        else
            txt = '  Q';
        end
        text(jpos(1),jpos(2),jpos(3),txt,'FontSize',20,'Color','k')
        
    end
    
    if S.n_sact>0
        
        CablePoints = S.CablePoints;
        XC_save     = CablePoints.XC_save;
        YC_save     = CablePoints.YC_save;
        ZC_save     = CablePoints.ZC_save;
        
        len_XC=length(XC_save);
        i1=1;
        for ii=1:len_XC
            if isnan(XC_save(ii))
                Pos = S.g_ini*[XC_save(i1:ii-1);YC_save(i1:ii-1);ZC_save(i1:ii-1);ones(size(XC_save(i1:ii-1)))];
                plot3(Pos(1,:),Pos(2,:),Pos(3,:),'Color','m');
                i1=ii+1;
            end
        end
    end
end


axis tight
xlim([-0.1*max(abs(xlim)),1.1*max(abs(xlim))])
ylim([-max(abs(ylim)),max(abs(ylim))])
zlim([-max(abs(zlim)),max(abs(zlim))])

drawnow
end