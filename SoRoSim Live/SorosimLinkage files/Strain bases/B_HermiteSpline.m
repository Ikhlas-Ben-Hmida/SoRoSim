function B = B_HermiteSpline(X,Bdof,Bodr)

% X varies from 0 to 1
% Bdof tells which deformation modes are enabled
% eg. Bdof = [1 1 1 0 0 0]' implies all 3 rotational modes are enabled
nele = max(Bodr); % Gives the number of elements
dof = sum(Bdof*(2*nele+2));
B   = zeros(6,dof);

k = 1;
for i=1:6
    w = 1/nele;
    a = 0;
        
    for j=1:Bdof(i)*nele
        b = a+w;
        if X>=a&&X<=b
            Xc = (X-a)/(b-a);
            B(i,k)   = 1-3*Xc^2+2*Xc^3;
            B(i,k+1) = Xc-2*Xc^2+Xc^3;
            B(i,k+2) = 3*Xc^2-2*Xc^3;
            B(i,k+3) = -Xc^2+Xc^3;
        end
        k = k+2;
        a = a+w;
    end
        
    k = k+Bdof(i);
end

end