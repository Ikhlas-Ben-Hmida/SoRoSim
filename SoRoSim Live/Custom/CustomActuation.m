function Fact = CustomActuation(Tr,q,g,J,t,qd,eta,Jdot)

%Tr: SorosimTree element, 
%q and qd: joint coordinates and their time derivatives, 
%g, J, Jd, and eta: transformation matrix, Jacobian, time derivative of jacobian, and screw velocity at every significant point of the linkage
%t:  time

%Fact should be 6*n column vector where n is the total number of gaussian
%points of all soft divisions including the start and end point of a division (nGauss).
%(Example: linkage with 2 soft links and 1 rigid link (n=nGauss1+nGauss2)
%Fact should be arranged according to the order of precedence of Link
%numbers

% Significant points: 1 for every joint, 1 at the center of the rigid link, 1
% at the start and end of every soft link and 1 for each Gaussian points

% J   = S.Jacobian(q);         %geometric jacobian of the linkage calculated at every significant points
% g   = S.FwdKinematics(q);    %transformation matrix of the linkage calculated at every significant points
% eta = S.ScrewVelocity(q,qd); %Screwvelocity of the linkage calculated at every significant points
% J   = S.Jacobiandot(q,qd);   %time derivative of geometric jacobian of the linkage calculated at every significant points


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CUSTOM INTERNAL ACTUATION, Tentacle actuation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=0;
for i=1:Tr.N
    for j=1:Tr.VLinks(i).npie-1
        n=n+Tr.CVTwists{i}(j+1).nip;
    end
end

%%

tf = 2;

Xs = Tr.CVTwists{1}(2).Xs;
Vect_M =zeros(1,length(Xs));
Vect_My =zeros(1,length(Xs));
Vect_Mz =zeros(1,length(Xs));

k = 1.45;
k = 1.36;
L0 = 0.20;

% if (t<=tf/2)
%     epsilon = 44.3+9.3*t*2;
% else
%     epsilon = 53.4-6*(t-1/2)*2;
% end

epsilon = 53.4;

epsilon = 53.4-10*t;

% t = 0;
E = 1e4;
% rbase = 0.005;
rbase = 0.01;
rtip = 0.001;
Xe = 0.5+0.4*t/tf;
% Xe = 0.6+0.3*t/tf;
Fact = [];

for ii = 1:length(Xs) 
    
    r(ii) = rbase - Xs(ii)*(rbase - rtip);
%     Jy(ii) = (pi/4)*(rbase^4);
    Jy(ii) = (pi/4)*(r(ii)^4);
    M(ii) = E*Jy(ii)*k*epsilon./(1+(epsilon*L0*(Xs(ii)-Xe))^2);
    strain(ii) = k*epsilon./(1+(epsilon*L0*(Xs(ii)-Xe))^2);
    Vect_M(ii) = M(ii);
    Fact = [Fact;[0 -M(ii) 0 0 0 0]'];
end

% Strain = k_force*Vect_M/(E*JJ);

% plot(Xs,Vect_M);
% a = 2;
% plot(Gauss,My,Gauss,Mz);legend('M_{y}','M_{z}')

% Fact = zeros(6*n,1); %change here

%%
end
