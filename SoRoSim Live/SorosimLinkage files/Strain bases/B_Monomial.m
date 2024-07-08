function B = B_Monomial(X,Bdof,Bodr)
dof = sum(Bdof.*(Bodr+1));
B   = zeros(6,dof);

k = 1;
for i=1:6
    for j=1:Bdof(i)*(Bodr(i)+1)
        B(i,k) = X^(j-1);
        k = k+1;
    end
end

end