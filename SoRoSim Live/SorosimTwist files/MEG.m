%This function computes Ms, Es, and Gs for a given Xs for the jth division
%of a link (15.12.2022)

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