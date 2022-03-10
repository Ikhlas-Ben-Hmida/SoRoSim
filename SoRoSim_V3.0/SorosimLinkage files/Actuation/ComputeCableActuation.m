%Function to compute cable actuation when cable is fully inside (14.06.2021)
function Bq = ComputeCableActuation(Tr,dc,dcp,Sdiv,Ediv,q)

Bq        = zeros(Tr.ndof,1);
dof_start = 1;

for i=1:Tr.N

    dof_start = dof_start+Tr.CVTwists{i}(1).dof;
    
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
       
        dof_here  = Tr.CVTwists{i}(j+1).dof;
        if isempty(dc{i}{j})
            dof_start = dof_start+dof_here;
            continue;
        end
            
        if j>=Sdiv{i}&&j<=Ediv{i}
            
            B         = Tr.CVTwists{i}(j+1).B;

            dcd  = dc{i}{j};
            dcpd = dcp{i}{j};
            
            q_here    = q(dof_start:dof_start+dof_here-1);
            xi_star   = Tr.CVTwists{i}(j+1).xi_star;
            ld        = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
            Ws        = Tr.VLinks(Tr.LinkIndex(i)).Ws{j};
            nGauss    = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};

            %scaling of quantities
            Lscale       = ld;
            ld           = 1;

            for k=2:nGauss-1

                dc_here     = dcd(:,k);
                dcp_here    = dcpd(:,k);
                B_here      = B((k-1)*6+1:k*6,:);
                xi_starhere = xi_star(6*(k-1)+1:6*k,1);
                xi_here     = B_here*q_here+xi_starhere;
                xihat_here  = dinamico_hat(xi_here);
                tang        = xihat_here*[dc_here;1]+[dcp_here;0];
                tang        = tang(1:3)/norm(tang(1:3)); %check with Federico
                Btau        = [dinamico_tilde(dc_here)*tang;tang];
                Btau(1:3)   = Btau(1:3)/Lscale;
                Bq(dof_start:dof_start+dof_here-1)  = Bq(dof_start:dof_start+dof_here-1)+ld*Ws(k)*B_here'*Btau*Lscale; %scaled back for addition
            end
            
        end
        dof_start = dof_start+dof_here;
    end 
end

end

