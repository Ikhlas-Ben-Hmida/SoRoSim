%Function to compute point force in Q space
%Last modified by Anup Teejo Mathew 02.03.2022
function Fp=ComputePointForce(Tr,J,g,t)

np = Tr.np;
Fp = zeros(Tr.ndof,1);
N  = Tr.N;

for ip=1:np
    done = false;
    Fp_loc = Tr.Fp_loc{ip};
    Fp_vec = Tr.Fp_vec{ip}(t);

    i_sig  = 1;

    for i=1:N

        i_sig = i_sig+1; %joint

        if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
            if i==Fp_loc(1)

                if ~Tr.FollowerForce{ip}
                    g_here        = g((i_sig-1)*4+1:i_sig*4,:);
                    g_here(1:3,4) = zeros(3,1);
                    Ad_g_here_inv = dinamico_Adjoint(ginv(g_here));
                    Fp_vec        = Ad_g_here_inv*Fp_vec;
                end

                Fp = Fp+J((i_sig-1)*6+1:i_sig*6,:)'*Fp_vec;
                break;
            end
            i_sig = i_sig+1;
        end

        for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
            Xs = Tr.CVTwists{i}(j+1).Xs;
            if i==Fp_loc(1)&&j==Fp_loc(2)
                Lscale      = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
                Fp_vec(1:3) = Fp_vec(1:3)/Lscale^2;
                Fp_vec(4:6) = Fp_vec(4:6)/Lscale;
                for ii=1:Tr.CVTwists{i}(j+1).nip
                    if Xs(ii)==Fp_loc(3)
                        if ~Tr.FollowerForce{ip}
                            g_here        = g((i_sig-1)*4+1:i_sig*4,:);
                            g_here(1:3,4) = zeros(3,1);
                            Ad_g_here_inv = dinamico_Adjoint(ginv(g_here));
                            Fp_vec        = Ad_g_here_inv*Fp_vec;
                        end
                        %scaling of quantities     
                        Fp = Fp+J((i_sig-1)*6+1:i_sig*6,:)'*Fp_vec*Lscale^2; %scaling back
                        done = true;
                        break;
                    end
                    i_sig = i_sig+1;
                end
                if done
                    break;
                end
            else
                i_sig = i_sig+Tr.CVTwists{i}(j+1).nip;
            end
        end
        if done
            break;
        end

    end

end

end

