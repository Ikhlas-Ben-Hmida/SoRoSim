%Generates a rotational symbol used by plotq0
%Last modified by Anup Teejo Mathew 03.02.2022
function Pos = RotationSymbol

z = 0.8;
y = 0.2;
z = [z;1];
y = [y;0];
z = [z;1.2];
y = [y;0.2];
z = [z;1];
y = [y;0];

theta = linspace(0,pi,18);
xarc  = cos(theta);
yarc  = sin(theta);

z = [z;xarc'];
y = [y;yarc'];
z = [z;-0.8];
y = [y;0.2];
z = [z;-1];
y = [y;0];
z = [z;-1.2];
y = [y;0.2];

Pos = [zeros(length(z),1) y z];

end

