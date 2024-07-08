function B = B_Chebychev(X,Bdof,Bodr)

% X varies from 0 to 1
% Bdof tells which deformation modes are enabled
% eg. Bdof = [1 1 1 0 0 0]' implies all 3 rotational modes are enabled
% Bodr defines the order of each strains
% eg. Bodr = [0 1 1 0 0 0]' implies linear bending strains while other are
% constant strains

dof = sum(Bdof.*(Bodr+1));
B   = zeros(6,dof);

k = 1;
X = 2*X-1; %transform from [0,1] to [-1 1]
for i=1:6
    T0 = 1;
    T1 = X;
    for j=1:Bdof(i)*(Bodr(i)+1)
        if j==1
            B(i,k) = T0;
        end
        if j==2
           B(i,k) = T1; 
        end
        if j>=3
            T1t = T1;
            T1 = 2*X*T1-T0;
            T0 = T1t;
            B(i,k) = T1;
        end
        k = k+1;
    end
end

end