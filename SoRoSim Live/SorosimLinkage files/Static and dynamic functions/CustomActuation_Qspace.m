%Function to convert user defined actuation force (Fact) on soft links to actuation force in Q space
%Last modified by Anup Teejo Mathew - 23/06/2021
function TauC = CustomActuation_Qspace(Tr,Fact)

ndof = Tr.ndof;
N    = Tr.N;

TauC    = zeros(ndof,1); %change here
i_sig_s = 1;

dof_start=1;
for i=1:N
    
    VTwists   = Tr.CVTwists{i};
    dof_start = dof_start+VTwists(1).dof;
    
    for j=1:Tr.VLinks(i).npie-1
        dof_here     = VTwists(j+1).dof;
        Lscale  = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Ws  = Tr.CVTwists{i}(j+1).Ws;
        nip = Tr.CVTwists{i}(j+1).nip;
        for ii=1:nip
            if Ws(ii)>0
                Fact_here      = Fact((i_sig_s-1)*6+1:i_sig_s*6);
                Fact_here(1:3) = Fact_here(1:3)/Lscale^2; %Nm
                Fact_here(4:6) = Fact_here(4:6)/Lscale;   %N

                B_here = VTwists(j+1).B((ii-1)*6+1:ii*6,:);

                TauC(dof_start:dof_start+dof_here-1) = TauC(dof_start:dof_start+dof_here-1)+Ws(ii)*B_here'*Fact_here*Lscale^2; %scaling back
            end
            i_sig_s = i_sig_s+1;
        end
        dof_start = dof_start+dof_here;
    end
end

end

