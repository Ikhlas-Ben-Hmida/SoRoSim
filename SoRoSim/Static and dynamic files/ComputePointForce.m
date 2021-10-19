%Function to compute point force in Q space
%Last modified by Anup Teejo Mathew - 24/05/2021

function Fp=ComputePointForce(S,J,t)

np     = S.np;
Fp     = zeros(S.ndof,1);
N      = S.N;

for ii=1:np
    
    Fp_loc = S.Fp_loc{ii};
    Fp_vec = S.Fp_vec{ii};
    Fp_vec = Fp_vec(t);
    i_sig  = 1;
    f      = 1;
    
    for i=1:N

        i_sig=i_sig+1;

        if f==Fp_loc
            Fp = Fp+J((i_sig-1)*6+1:i_sig*6,:)'*Fp_vec;
        end

        if S.VLinks(i).linktype=='r'
            i_sig = i_sig+1;
        end
        f=f+1;

        if f>Fp_loc
            break;
        end

        for j=1:S.VLinks(i).npie-1

            %scaling of quantities
            Lscale      = S.VLinks(i).lp{j};
            Fp_vec(1:3) = Fp_vec(1:3)/Lscale^2;
            Fp_vec(4:6) = Fp_vec(4:6)/Lscale;

            i_sig = i_sig+S.VLinks(i).nGauss{j}-1;
            
            if f==Fp_loc
                Fp=Fp+J((i_sig-1)*6+1:i_sig*6,:)'*Fp_vec*Lscale^2; %scaling back 
            end
            
            i_sig = i_sig+1;
            f     = f+1;

            if f>Fp_loc
                break;
            end
            
        end

    end

end

end

