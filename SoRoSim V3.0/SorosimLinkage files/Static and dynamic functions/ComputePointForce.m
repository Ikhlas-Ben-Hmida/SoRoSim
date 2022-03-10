%Function to compute point force in Q space
%Last modified by Anup Teejo Mathew 02.03.2022
function Fp=ComputePointForce(Tr,J,t,g)

np = Tr.np;
Fp = zeros(Tr.ndof,1);
N  = Tr.N;

for ii=1:np

    Fp_loc = Tr.Fp_loc{ii};
    Fp_vec = Tr.Fp_vec{ii};
    Fp_vec = Fp_vec(t);
    i_sig  = 1;

    for i=1:N

        i_sig = i_sig+1; %joint

        if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
            if i==Fp_loc(1)

                if ~Tr.FollowerForce{ii}
                    g_here        = g((i_sig-1)*4+1:i_sig*4,:);
                    g_here(1:3,4) = zeros(3,1);
                    Adg_here      = dinamico_Adjoint(ginv(g_here));
                    Fp_vec        = Adg_here*Fp_vec;
                end

                Fp = Fp+J((i_sig-1)*6+1:i_sig*6,:)'*Fp_vec;
                break;
            end
            i_sig = i_sig+1;
        end

        for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1

            i_sig = i_sig+Tr.VLinks(Tr.LinkIndex(i)).nGauss{j}-1;

            if i==Fp_loc(1)&&j==Fp_loc(2)
                Lscale      = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
                Fp_vec(1:3) = Fp_vec(1:3)/Lscale^2;
                Fp_vec(4:6) = Fp_vec(4:6)/Lscale;

                if ~Tr.FollowerForce{ii}
                    g_here        = g((i_sig-1)*4+1:i_sig*4,:);
                    g_here(1:3,4) = zeros(3,1);
                    Adg_here      = dinamico_Adjoint(ginv(g_here));
                    Fp_vec        = Adg_here*Fp_vec;
                end

                %scaling of quantities     
                Fp = Fp+J((i_sig-1)*6+1:i_sig*6,:)'*Fp_vec*Lscale^2; %scaling back
                break;
            end

            i_sig = i_sig+1;

        end

    end

end

end

