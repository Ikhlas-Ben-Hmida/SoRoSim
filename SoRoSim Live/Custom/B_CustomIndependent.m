function B = B_CustomIndependent(X,varargin)
%To add a different custom base, make a file similar to this and change the file name of function handle in the corresponding Twist
%column of B should be linearly independent with others
%Change SorosimTwist.Type = 'Custom Independent'

% X varies from 0 to 1

%Two FEM Quad Elements with constant torsion and constant elongation

nele = 1;

dof = 18;

B   = zeros(6,dof);

k = 1;

Bdof = [1 1 1 1 1 1]'; %FEM dof

for i=1:6
    w = 1/nele;
    a = 0;

    for j=1:Bdof(i)*nele
        b = a+w;
        if X>=a&&X<=b
            Xc = (X-a)/(b-a);
            B(i,k)   =  9/2*(Xc-1/2)*(Xc-5/6);
            B(i,k+1) = -9*(Xc-1/6)*(Xc-5/6);
            B(i,k+2) =  9/2*(Xc-1/6)*(Xc-1/2);
        end
        k = k+2;
        a = a+w;
    end

    k = k+Bdof(i);
end

end