function u = CustomActuatorStrength(S,q,g,J,t,qd,eta,Jd)

%S: Linkage element, 
%q and qd: joint coordinates and their time derivatives, 
%g, J, Jd, and eta: transformation matrix, Jacobian, time derivative of jacobian, and screw velocity at every significant point of the linkage
%t:  time

% Significant points: 1 for every joint, 1 at the center of the rigid link, 1
% at the start and end of every soft link and 1 for each Gauss quadrature points

% J   = S.Jacobian(q);         %geometric jacobian of the linkage calculated at every significant points
% g   = S.FwdKinematics(q);    %transformation matrix of the linkage calculated at every significant points
% eta = S.ScrewVelocity(q,qd); %Screwvelocity of the linkage calculated at every significant points
% J   = S.Jacobiandot(q,qd);   %time derivative of geometric jacobian of the linkage calculated at every significant points

nact = S.nact;
u = zeros(nact,1); %change here

end

