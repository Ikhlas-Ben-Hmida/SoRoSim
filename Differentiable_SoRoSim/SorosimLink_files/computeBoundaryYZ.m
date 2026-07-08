%Function to compute boundary points of a cross section
%Last modified by Anup Teejo Mathew 28.11.2024
%Modified to support the 'am_isupport' cross section. To keep plotq's
%patch stitching valid, the boundary is a single closed loop of n_r points
%(same convention as 'C' and 'E'). The AM I-Support outline is rendered as
%the body envelope with three shallow lobes at the chamber angles theta.
%A companion function computeChamberYZ (below) returns the individual
%annular chamber outlines for line-based plotting (used by plotq0).
function [y,z] = computeBoundaryYZ(Link,X,varargin) % X varies from 0 to 1, varargin is division number (only for soft link)

n_r = Link.n_r;

if Link.CS=='C'

    if Link.linktype == 'r'
        r_fn  = Link.r;
    else
        j = varargin{1}; %division number
        r_fn  = Link.r{j};
    end
    r     = r_fn(X);
    theta = linspace(0,2*pi,n_r);
    y     = r*sin(theta);
    z     = r*cos(theta);

elseif Link.CS=='R'

    if Link.linktype == 'r'
        h_fn  = Link.h;
        w_fn  = Link.w;
    else
        j = varargin{1}; %division number
        h_fn  = Link.h{j};
        w_fn  = Link.w{j};
    end
    h     = h_fn(X);
    w     = w_fn(X);
    y     = [h/2 -h/2 -h/2 h/2 h/2];
    z     = [w/2 w/2 -w/2 -w/2 w/2];

elseif Link.CS=='E'

    if Link.linktype == 'r'
        a_fn  = Link.a;
        b_fn  = Link.b;
    else
        j = varargin{1}; %division number
        a_fn  = Link.a{j};
        b_fn  = Link.b{j};
    end
    a     = a_fn(X);
    b     = b_fn(X);
    theta = linspace(0,2*pi,n_r);
    y     = a*sin(theta);
    z     = b*cos(theta);

elseif strcmp(Link.CS,'am_isupport')

    %Link.r stores the body envelope radius (function of X1), set at
    %construction. delta and theta describe the chamber layout.
    j     = varargin{1}; %division number (am_isupport is soft only)
    r_fn  = Link.r{j};
    r_env = r_fn(X);

    ang   = linspace(0,2*pi,n_r);

    %Optional three-lobe modulation so the outline suggests the chambers.
    %Purely cosmetic: bumps the envelope outward near each chamber angle.
    %Set lobe_amp = 0 to fall back to a plain circle of radius r_env.
    lobe_amp = 0.08*r_env; %8% of the envelope radius
    ro_fn    = Link.ro{j};
    ro_here  = ro_fn(X);
    th_c     = Link.theta; %chamber angles [rad], measured from local y (d1 in paper)

    bump = zeros(size(ang));
    for k=1:numel(th_c)
        %angular distance to chamber k (wrapped to [-pi,pi])
        dth  = atan2(sin(ang-th_c(k)),cos(ang-th_c(k)));
        %narrow raised-cosine bump centred on the chamber
        sigma = max(ro_here/r_env, 0.15); %angular width tied to chamber size
        bump  = bump + exp(-(dth.^2)/(2*sigma^2));
    end
    r_loop = r_env + lobe_amp*bump;

    %theta measured from d1 (local y): y along d1, z along d2
    y = r_loop.*cos(ang);
    z = r_loop.*sin(ang);

end

end