%Function to convert user defined actuation force (Fact) on soft links to actuation force in Q space
%Last modified by Anup Teejo Mathew - 23/06/2021
function TauC = CustomActuation_Qspace(Tr,Fact)

ndof = Tr.ndof;
N    = Tr.N;

TauC    = zeros(ndof,1); %change here
i_sig_s = 1;

dof_start=1;
for i=1:N
    dof_start = dof_start+Tr.CVtwists{i}(1).dof;
    for j=1:Tr.VLinks(i).npie-1
        B       = Tr.CVtwists{i}(j+1).B;
        dof     = Tr.CVtwists{i}(j+1).dof;
        i_sig_s = i_sig_s+1;
        Lscale  = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Ws      = Tr.VLinks(Tr.LinkIndex(i)).Ws{j};
        for k=2:Tr.VLinks(Tr.LinkIndex(i)).nGauss{j}-1
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

