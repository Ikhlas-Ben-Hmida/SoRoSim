%Function to compute cable actuation when cable is fully inside (14.06.2021)
function Bq = ComputeCableActuation(S,dc,dcp,Sdiv,Ediv,q)

Bq        = zeros(S.ndof,1);
f         = 1;
dof_start = 1;

for i=1:S.N

    dof_start = dof_start+S.Vtwists(f).dof;
    f         = f+1;
    dci       = dc{i};
    dcpi      = dcp{i};
    Sdivi     = Sdiv{i};
    Edivi     = Ediv{i};
    for j=1:S.VLinks(i).npie-1

        B         = S.Vtwists(f).B;
        dof_here  = S.Vtwists(f).dof;
            
        if j>=Sdivi&&j<=Edivi

            dcd  = dci{j};
            dcpd = dcpi{j};
            
            q_here    = q(dof_start:dof_start+dof_here-1);
            xi_star   = S.Vtwists(f).xi_star;
            lpf       = S.VLinks(i).lp{j};
            Ws        = S.VLinks(i).Ws{j};
            nGauss    = S.VLinks(i).nGauss{j};

            %scaling of quantities
            Lscale        = lpf;
            lpf           = 1;

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
                Bq(dof_start:dof_start+dof_here-1)  = Bq(dof_start:dof_start+dof_here-1)+lpf*Ws(k)*B_here'*Btau*Lscale; %scaled back for addition
            end
            
        end
        dof_start = dof_start+dof_here;
        f = f+1;
    end 
end

end

