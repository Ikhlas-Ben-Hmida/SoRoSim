%This function updates dependent properties when a property is changed
%(20.05.2021)

function [varargout]= PropUpdate(Link)


%% Rigid Link
if Link.linktype=='r'

    if Link.CS == 'R'
        %unaffected properties:
        h_fn = Link.h;
        w_fn = Link.w;
        L    = Link.L;
        Rho  = Link.Rho;

        %Updating properties:
        [Xs,Ws,nGauss] = GaussQuadrature(10);
        V      = 0;
        M_area = 0;

        for ii=2:nGauss-1
            h_here = h_fn(Xs(ii));
            w_here = w_fn(Xs(ii));
            A_here = h_here*w_here;
            V      = V + Ws(ii)*L*A_here;
            M_area = M_area+ Ws(ii)*L*A_here*Xs(ii)*L;
        end

        cx   = M_area/V;
        mass = V*Rho;

        Iy  = 0;
        Iz  = 0;
        for ii=2:nGauss-1
            h_here = h_fn(Xs(ii));
            w_here = w_fn(Xs(ii));
            Iy     = Iy+Ws(ii)*L*((Xs(ii)*L-cx)^2*h_here*w_here+h_here*w_here^3/12);
            Iz     = Iz+Ws(ii)*L*((Xs(ii)*L-cx)^2*h_here*w_here+h_here^3*w_here/12);
        end
        Ix = Iy+Iz;
    else
        %unaffected properties:
        r_fn = Link.r;
        Rho  = Link.Rho;
        L    = Link.L;

        %Updating properties:
        [Xs,Ws,nGauss] = GaussQuadrature(10);
        V      = 0;
        M_area = 0;

        for ii=2:nGauss-1
            r_here = r_fn(Xs(ii));
            A_here = pi*r_here^2;
            V      = V + Ws(ii)*L*A_here;
            M_area = M_area+ Ws(ii)*L*A_here*Xs(ii)*L;
        end

        cx   = M_area/V;
        mass = V*Rho;

        Iy  = 0;
        for ii=2:nGauss-1
            r_here = r_fn(Xs(ii));
            Iy     = Iy+Ws(ii)*L*((pi/4)*r_here^4+pi*r_here^2*(Xs(ii)*L-cx)^2);
        end
        Iz = Iy;
        Ix = Iy+Iz;
    end
    
    Link.cx = cx;
    
    Ms = double([Ix 0 0 0 0 0;0 Iy 0 0 0 0;0 0 Iz 0 0 0;...
        0 0 0 mass 0 0;0 0 0 0 mass 0;0 0 0 0 0 mass]);

    varargout{1} = Ms; %Ms
    
else
    %% Soft Link
    ndiv         = Link.npie-1;
    Ms           = cell(1,ndiv); Es = cell(1,ndiv); Gs = cell(1,ndiv);

    for i=1:(Link.npie)-1

        if Link.CS=='R'
            h_fn     = Link.h{i};
            w_fn     = Link.w{i};
            Rho      = Link.Rho;
            G        = Link.G;
            E        = Link.E;
            Eta      = Link.Eta;
            Xs_p     = Link.Xs{i};
            nGauss_p = Link.nGauss{i};
            l_p      = Link.lp{i};

            %updating:
            h_nGauss = zeros(nGauss_p,1);
            w_nGauss = zeros(nGauss_p,1);
            for ii=1:nGauss_p
                h_nGauss(ii) = h_fn(Xs_p(ii));
                w_nGauss(ii) = w_fn(Xs_p(ii));
            end

            Iy_p = (1/12)*h_nGauss.*(w_nGauss.^3);
            Iz_p = (1/12)*(h_nGauss.^3).*w_nGauss;
            Ix_p = Iy_p+Iz_p;
            A_p  = h_nGauss.*w_nGauss;

        else %Circular
            r_fn     = Link.r{i};
            Rho      = Link.Rho;
            G        = Link.G;
            E        = Link.E;
            Eta      = Link.Eta;
            Xs_p     = Link.Xs{i};
            nGauss_p = Link.nGauss{i};
            l_p      = Link.lp{i};

            %updating:
            r_nGauss = zeros(nGauss_p,1);
            for ii=1:nGauss_p
                r_nGauss(ii) = r_fn(Xs_p(ii));
            end
            Iy_p = (pi/4)*r_nGauss.^4;
            Iz_p = Iy_p;
            Ix_p = Iy_p+Iz_p;
            A_p  = pi*r_nGauss.^2;

        end

        Ms_p = zeros(6*nGauss_p,6); %inertia
        Es_p = zeros(6*nGauss_p,6); %stiffness
        Gs_p = zeros(6*nGauss_p,6); %damping

        for ii=1:nGauss_p
            Ms_p((ii-1)*6+1:ii*6,:) = Rho*diag([Ix_p(ii),Iy_p(ii),Iz_p(ii),A_p(ii),A_p(ii),A_p(ii)]);
            Es_p((ii-1)*6+1:ii*6,:) = diag([G*Ix_p(ii),E*Iy_p(ii),E*Iz_p(ii),E*A_p(ii),G*A_p(ii),G*A_p(ii)]);
            Gs_p((ii-1)*6+1:ii*6,:) = Eta*diag([Ix_p(ii),3*Iy_p(ii),3*Iz_p(ii),3*A_p(ii),A_p(ii),A_p(ii)]);
        end

        Ms{i} = Ms_p;
        Es{i} = Es_p;
        Gs{i} = Gs_p;
    end

    varargout{1} = Ms; %Ms
    varargout{2} = Es; %Es
    varargout{3} = Gs; %Gs

    
end

%eof