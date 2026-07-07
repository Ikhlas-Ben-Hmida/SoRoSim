%This function computes Ms, Es, and Gs for a given Xs for the jth division
%of a link (15.12.2022)
%Modified to support the 'am_isupport' cross section: N annular pneumatic
%chambers at radial distance delta from the centroid, disposed at angular
%positions theta (Alessi, Falotico, Lucantonio, IEEE Access 2023,
%doi: 10.1109/ACCESS.2023.3266282, Eqs. (33)-(34))

function [Ms,Es,Gs]= MEG(Link,j,Xs)

np = length(Xs);

%j is division number
if Link.CS=='R'

    h_fn     = Link.h{j};
    w_fn     = Link.w{j};

    %updating:
    h_nGauss = zeros(np,1);
    w_nGauss = zeros(np,1);

    for ii=1:np
        h_nGauss(ii) = h_fn(Xs(ii));
        w_nGauss(ii) = w_fn(Xs(ii));
    end

    Iy_p = (1/12)*h_nGauss.*(w_nGauss.^3);
    Iz_p = (1/12)*(h_nGauss.^3).*w_nGauss;
    Ix_p = Iy_p+Iz_p;
    A_p  = h_nGauss.*w_nGauss;

elseif Link.CS=='C' %Circular

    r_fn     = Link.r{j};

    %updating:
    r_nGauss = zeros(np,1);

    for ii=1:np
        r_nGauss(ii) = r_fn(Xs(ii));
    end

    Iy_p = (pi/4)*r_nGauss.^4;
    Iz_p = Iy_p;
    Ix_p = Iy_p+Iz_p;
    A_p  = pi*r_nGauss.^2;

elseif Link.CS=='E'

    a_fn     = Link.a{j};
    b_fn     = Link.b{j};

    %updating:
    a_nGauss = zeros(np,1);
    b_nGauss = zeros(np,1);

    for ii=1:np
        a_nGauss(ii) = a_fn(Xs(ii));
        b_nGauss(ii) = b_fn(Xs(ii));
    end

    Iy_p = (pi/4)*a_nGauss.*(b_nGauss.^3);
    Iz_p = (pi/4)*(a_nGauss.^3).*b_nGauss;
    Ix_p = Iy_p+Iz_p;
    A_p  = pi*a_nGauss.*b_nGauss;

elseif strcmp(Link.CS,'am_isupport') %AM I-Support: N annular pneumatic chambers

    ro_fn = Link.ro{j};       %outer radius of one chamber, function of X1
    ri_fn = Link.ri{j};       %inner radius of one chamber, function of X1
    delta = Link.delta;       %radial distance of chamber centres [m]
    theta = Link.theta;       %(1xN) angular positions of chambers [rad]
    N_ch  = length(theta);    %number of chambers

    %updating:
    ro_nGauss = zeros(np,1);
    ri_nGauss = zeros(np,1);

    for ii=1:np
        ro_nGauss(ii) = ro_fn(Xs(ii));
        ri_nGauss(ii) = ri_fn(Xs(ii));
    end

    a_p  = pi*(ro_nGauss.^2-ri_nGauss.^2);      %annulus area of a single chamber
    A_p  = N_ch*a_p;                            %effective cross-section area: A = N*a (disks ignored)
    I0_p = (pi/4)*(ro_nGauss.^4-ri_nGauss.^4);  %second moment of area of one annulus about its own centroid

    %Second moments of area about the local transverse axes (parallel axis
    %theorem), Eqs. (33)-(34) of the paper:
    %I1 = sum_j [ pi/4 (ro^4-ri^4) + a (delta*sin(theta_j))^2 ]
    %I2 = sum_j [ pi/4 (ro^4-ri^4) + a (delta*cos(theta_j))^2 ]
    Iy_p = zeros(np,1);
    Iz_p = zeros(np,1);
    for jj=1:N_ch
        Iy_p = Iy_p+I0_p+a_p*(delta*sin(theta(jj)))^2;
        Iz_p = Iz_p+I0_p+a_p*(delta*cos(theta(jj)))^2;
    end
    Ix_p = Iy_p+Iz_p; %polar moment of area: I3 = I1+I2

end

Ms = zeros(6*np,6); %inertia
Es = zeros(6*np,6); %stiffness
Gs = zeros(6*np,6); %damping

Rho  = Link.Rho;
G    = Link.G;
E    = Link.E;
Eta  = Link.Eta;

for ii=1:np
    Ms((ii-1)*6+1:ii*6,:) = Rho*diag([Ix_p(ii),Iy_p(ii),Iz_p(ii),A_p(ii),A_p(ii),A_p(ii)]);
    Es((ii-1)*6+1:ii*6,:) = diag([G*Ix_p(ii),E*Iy_p(ii),E*Iz_p(ii),E*A_p(ii),G*A_p(ii),G*A_p(ii)]);
    Gs((ii-1)*6+1:ii*6,:) = Eta*diag([Ix_p(ii),3*Iy_p(ii),3*Iz_p(ii),3*A_p(ii),A_p(ii),A_p(ii)]);
end

%eof