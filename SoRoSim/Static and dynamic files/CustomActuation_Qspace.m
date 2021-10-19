%Function to convert user defined actuation force (Fact) on soft links to actuation force in Q space
%Last modified by Anup Teejo Mathew - 23/06/2021
function TauC = CustomActuation_Qspace(S,Fact)

ndof = S.ndof;
N = S.N;

TauC = zeros(ndof,1); %change here
i_sig_s = 1;

f=1;
dof_start=1;
for i=1:N
    f=f+1;
    for j=1:S.VLinks(i).npie-1
        B       = S.Vtwists(f).B;
        dof     = S.Vtwists(f).dof;
        i_sig_s = i_sig_s+1;
        Lscale  = S.VLinks(i).lp{j};
        Ws      = S.VLinks(i).Ws{j};
        for k=2:S.VLinks(i).nGauss{j}-1
            Fact_here      = Fact((i_sig_s-1)*6+1:i_sig_s*6);
            Fact_here(1:3) = Fact_here(1:3)/Lscale^2; %Nm
            Fact_here(4:6) = Fact_here(4:6)/Lscale;   %N
            B_here         = B((k-1)*6+1:k*6,:);

            TauC(dof_start:dof_start+dof-1) = TauC(dof_start:dof_start+dof-1)+Ws(k)*B_here'*Fact_here*Lscale^2; %scaling back
        end
        dof_start = dof_start+dof;
    end
end

end

