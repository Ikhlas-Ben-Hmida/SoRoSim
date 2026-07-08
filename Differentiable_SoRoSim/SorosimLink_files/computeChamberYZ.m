%Outlines of the individual annular pneumatic chambers for the
%'am_isupport' cross section. Returns cell arrays Yc, Zc where each cell is
%one chamber outline (outer circle, a NaN break, then the inner circle) so a
%single plot3 call draws the annulus without a connecting chord.
%Intended for line plotting (plot3), NOT for patch stitching.
%Chamber centre k sits at (delta*cos(theta_k), delta*sin(theta_k)) in the
%(d1,d2) = (y,z) plane, matching MEG.m and the paper (Eqs. 33-34, Fig. 2).
%Last modified 08.07.2025
function [Yc,Zc] = computeChamberYZ(Link,X,j)

if ~strcmp(Link.CS,'am_isupport')
    Yc = {}; Zc = {};
    return
end

n_r   = Link.n_r;
ang   = linspace(0,2*pi,n_r);

ro    = Link.ro{j}(X); %outer chamber radius at this station
ri    = Link.ri{j}(X); %inner chamber radius at this station
delta = Link.delta;    %radial distance of chamber centres
th_c  = Link.theta;    %chamber angles [rad] measured from d1 (y)

N_ch  = numel(th_c);
Yc    = cell(1,N_ch);
Zc    = cell(1,N_ch);

for k=1:N_ch
    yc = delta*cos(th_c(k)); %chamber centre (y along d1)
    zc = delta*sin(th_c(k)); %chamber centre (z along d2)

    y_out = yc + ro*cos(ang);
    z_out = zc + ro*sin(ang);
    y_in  = yc + ri*cos(ang);
    z_in  = zc + ri*sin(ang);

    Yc{k} = [y_out, NaN, y_in];
    Zc{k} = [z_out, NaN, z_in];
end

end