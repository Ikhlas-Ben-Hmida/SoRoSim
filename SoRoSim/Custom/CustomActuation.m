function Fact = CustomActuation(S,q,g,J,t,qd,eta,Jdot)


%S: Linkage element, 
%q and qd: joint coordinates and their time derivatives, 
%g, J, Jd, and eta: transformation matrix, Jacobian, time derivative of jacobian, and screw velocity at every significant point of the linkage
%t:  time

%Fact should be 6*n column vector where n is the total number of gaussian
%points of all soft divisions including the start and end point of a division (nGauss).
%(Example: linkage with 2 soft links and 1 rigid link (n=nGauss1+nGauss2)
%Fact should be arranged according to the order of precedence

% Significant points: 1 for every joint, 1 at the center of the rigid link, 1
% at the start and end of every soft link and 1 for each Gaussian points

% J   = S.Jacobian(q);         %geometric jacobian of the linkage calculated at every significant points
% g   = S.FwdKinematics(q);    %transformation matrix of the linkage calculated at every significant points
% eta = S.ScrewVelocity(q,qd); %Screwvelocity of the linkage calculated at every significant points
% J   = S.Jacobiandot(q,qd);   %time derivative of geometric jacobian of the linkage calculated at every significant points

n=0;
for i=1:S.N
    for j=1:S.VLinks(i).npie-1
        n=n+S.VLinks(i).nGauss(j);
    end
end
Fact = zeros(6*n,1); %change here

end

