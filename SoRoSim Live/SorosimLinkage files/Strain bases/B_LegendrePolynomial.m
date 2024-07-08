function B = B_LegendrePolynomial(X,Bdof,Bodr)

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
    P0 = 1;
    P1 = X;
    for j=1:Bdof(i)*(Bodr(i)+1)
        if j==1
            B(i,k) = P0;
        end
        if j==2
           B(i,k) = P1; 
        end
        if j>=3
            n=j-2;
            P1t = P1;
            P1 = ((2*n+1)*X*P1-n*P0)/(n+1);
            P0 = P1t;
            B(i,k) = P1;
        end
        k = k+1;
    end
end


end