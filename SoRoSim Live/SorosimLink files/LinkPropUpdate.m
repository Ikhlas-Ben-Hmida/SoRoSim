%This function updates dependent properties when a property is changed
%(15.12.2021)

function M = LinkPropUpdate(Link)

%% Rigid Link only

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

    if V~=0
        cx   = M_area/V;
    else
        if L==0
            cx = 0;
        else
            Vh      = 0;
            Mh_area = 0;
            Vw      = 0;
            Mw_area = 0;
            for ii=2:nGauss-1
                h_here = h_fn(Xs(ii));
                w_here = w_fn(Xs(ii));
                Vh      = Vh + Ws(ii)*L*h_here;
                Mh_area = Mh_area+ Ws(ii)*L*h_here*Xs(ii)*L;
                Vw      = Vw + Ws(ii)*L*w_here;
                Mw_area = Mw_area+ Ws(ii)*L*w_here*Xs(ii)*L;
            end
            if Vh==0
                cx   = Mw_area/Vw;
            else
                cx   = Mh_area/Vh;
            end
        end
    end

    mass = V*Rho;

    Ix = 0;
    Iy  = 0;
    Iz  = 0;
    for ii=2:nGauss-1
        h_here = h_fn(Xs(ii));
        w_here = w_fn(Xs(ii));
        Iy     = Iy+Ws(ii)*L*((Xs(ii)*L-cx)^2*h_here*w_here+h_here*w_here^3/12);
        Iz     = Iz+Ws(ii)*L*((Xs(ii)*L-cx)^2*h_here*w_here+h_here^3*w_here/12);
        Ix     = Iz+Ws(ii)*L*(h_here*w_here^3/12+h_here^3*w_here/12);
    end

elseif Link.CS=='C'
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

    if V~=0
        cx   = M_area/V;
    else
        if L==0
            cx = 0;
        else
            cx = L/2;
        end
    end
    mass = V*Rho;

    Ix = 0;
    Iy = 0;
    for ii=2:nGauss-1
        r_here = r_fn(Xs(ii));
        Iy     = Iy+Ws(ii)*L*((pi/4)*r_here^4+pi*r_here^2*(Xs(ii)*L-cx)^2);
        Ix     = Ix+Ws(ii)*L*((pi/2)*r_here^4);
    end
    Iz = Iy;

elseif Link.CS=='E'
    %unaffected properties:
    a_fn = Link.a;
    b_fn = Link.b;
    L    = Link.L;
    Rho  = Link.Rho;

    %Updating properties:
    [Xs,Ws,nGauss] = GaussQuadrature(10);
    V      = 0;
    M_area = 0;

    for ii=2:nGauss-1
        a_here = a_fn(Xs(ii));
        b_here = b_fn(Xs(ii));
        A_here = pi*a_here*b_here;
        V      = V + Ws(ii)*L*A_here;
        M_area = M_area+ Ws(ii)*L*A_here*Xs(ii)*L;
    end

    if V~=0
        cx   = M_area/V;
    else
        if L==0
            cx = 0;
        else
            Va      = 0;
            Ma_area = 0;
            Vb      = 0;
            Mb_area = 0;
            for ii=2:nGauss-1
                a_here = a_fn(Xs(ii));
                b_here = b_fn(Xs(ii));
                Va      = Va + Ws(ii)*L*a_here;
                Ma_area = Ma_area+ Ws(ii)*L*a_here*Xs(ii)*L;
                Vb      = Vb + Ws(ii)*L*b_here;
                Mb_area = Mb_area+ Ws(ii)*L*b_here*Xs(ii)*L;
            end
            if Va==0
                cx   = Mb_area/Vb;
            else
                cx   = Ma_area/Va;
            end
        end
    end
    mass = V*Rho;

    Ix = 0;
    Iy  = 0;
    Iz  = 0;
    for ii=2:nGauss-1
        a_here = a_fn(Xs(ii));
        b_here = b_fn(Xs(ii));
        Iy     = Iy+Ws(ii)*L*((Xs(ii)*L-cx)^2*pi*a_here*b_here+pi*a_here*b_here^3/4);
        Iz     = Iz+Ws(ii)*L*((Xs(ii)*L-cx)^2*pi*a_here*b_here+pi*a_here^3*b_here/4);
        Ix     = Ix+Ws(ii)*L*(pi*a_here*b_here^3/4+pi*a_here^3*b_here/4);
    end
end

Ix = Ix*Rho;
Iy = Iy*Rho;
Iz = Iz*Rho;

M = double([Ix 0 0 0 0 0;0 Iy 0 0 0 0;0 0 Iz 0 0 0;...
             0 0 0 mass 0 0;0 0 0 0 mass 0;0 0 0 0 0 mass]);

Link.gi = [eye(3),[cx;0;0];[0 0 0 1]];
Link.gf = [eye(3),[L-cx;0;0];[0 0 0 1]];

    
end

%eof