function B = B_IMQ(X,Bdof,Bodr)

% X varies from 0 to 1
% Bdof tells which deformation modes are enabled
% eg. Bdof = [1 1 1 0 0 0]' implies all 3 rotational modes are enabled
% Bodr defines the order of each strains
% eg. Bodr = [0 1 1 0 0 0]' implies linear bending strains while other are
% constant strains

dof = sum(Bdof.*(Bodr+1));
B   = zeros(6,dof);

k = 1;
for i=1:6
    if Bodr(i)==0&&Bdof(i)==1
        B(i,k) = 1;
        k = k+1;
    else
        w = 1/Bodr(i);
        a = 0;
        c = 2*sqrt(3)/w; % to have a value of 0.5 at w/2
        for j=1:Bdof(i)*(Bodr(i)+1)
            B(i,k) = 1/sqrt(1+(X-a)^2*c^2);
            k = k+1;
            a = a+w;
        end
    end
end

end