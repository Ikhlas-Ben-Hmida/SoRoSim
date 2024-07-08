function B = B_FEMLike(X,Bdof,Bodr,FEMOrder)

% X varies from 0 to 1
% Bdof tells which deformation modes are enabled
% eg. Bdof = [1 1 1 0 0 0]' implies all 3 rotational modes are enabled
% Linear, Quadratic and Cubic FEM basis
nele = max(Bodr); % Gives the number of elements

if strcmp(FEMOrder,'Linear')
    dof = sum(Bdof*(nele+1));
elseif strcmp(FEMOrder,'Quadratic')
    dof = sum(Bdof*(2*nele+1));
else % cubic
    dof = sum(Bdof*(3*nele+1));
end
B   = zeros(6,dof);

k = 1;
for i=1:6
    w = 1/nele;
    a = 0;
    if strcmp(FEMOrder,'Linear')

        for j=1:Bdof(i)*nele
            b = a+w;
            if X>=a&&(X<b||(X-b)^2<1e-12)
                Xc = (X-a)/(b-a);
                B(i,k) = 1-Xc;
                B(i,k+1) = Xc;
            end
            k = k+1;
            a = a+w;
        end
        
    elseif strcmp(FEMOrder,'Quadratic')

        for j=1:Bdof(i)*nele
            b = a+w;
            if X>=a&&X<=b
                Xc = (X-a)/(b-a);
                B(i,k)   =  2*(Xc-1/2)*(Xc-1);
                B(i,k+1) = -4*Xc*(Xc-1);
                B(i,k+2) =  2*Xc*(Xc-1/2);
            end
            k = k+2;
            a = a+w;
        end
        
    else %Cubic
        
        for j=1:Bdof(i)*nele
            b = a+w;
            if X>=a&&X<=b
                Xc = (X-a)/(b-a);
                B(i,k)   = -9/2*(Xc-1/3)*(Xc-2/3)*(Xc-1);
                B(i,k+1) =  27/2*Xc*(Xc-2/3)*(Xc-1);
                B(i,k+2) = -27/2*Xc*(Xc-1/3)*(Xc-1);
                B(i,k+3) =  9/2*Xc*(Xc-1/3)*(Xc-2/3);
            end
            k = k+3;
            a = a+w;
        end
        
    end
    k = k+Bdof(i);
end

end