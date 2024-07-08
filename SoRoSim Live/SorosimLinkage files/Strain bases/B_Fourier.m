function B = B_Fourier(X,Bdof,Bodr)
dof = sum(Bdof.*(2*Bodr+1));
B   = zeros(6,dof);

k = 1;
for i=1:6
    for j=1:Bdof(i)*(Bodr(i)+1)
        if j==1
            B(i,k) = 1;
            k = k+1;
        else
            B(i,k)   = cos(2*pi*(j-1)*X);
            B(i,k+1) = sin(2*pi*(j-1)*X);
            k = k+2;
        end
    end
end

end