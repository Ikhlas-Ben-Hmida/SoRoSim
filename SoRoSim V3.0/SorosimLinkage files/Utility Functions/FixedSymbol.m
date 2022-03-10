%Generates a fixed symbol used by plotq0 to denote fixed joints
%Last modified by Anup Teejo Mathew 03.02.2022

function Pos = FixedSymbol

y = 1;
z = 1;
y = [y;-1];
z = [z;1];
y = [y;-1];
z = [z;-1];
y = [y;1];
z = [z;-1];
y = [y;1];
z = [z;1];
y = [y;-1];
z = [z;-1];
y = [y;-1];
z = [z;1];
y = [y;1];
z = [z;-1];

Pos = [zeros(length(z),1) y z];

end

