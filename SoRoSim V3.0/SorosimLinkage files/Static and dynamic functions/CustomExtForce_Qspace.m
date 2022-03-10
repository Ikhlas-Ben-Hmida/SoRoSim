%Function to convert user defined external force (Fext) on each significant
%points to external force in Q space
%Last modified by Anup Teejo Mathew - 24/05/2021
function [F_custom,M_custom,C_custom]=CustomExtForce_Qspace(Tr,J,Fext,Jd,eta)


M_added  = Tr.M_added;
N        = Tr.N;
ndof     = Tr.ndof;

n_sig    = Tr.nsig;

if nargin==3
    Jd  = zeros(6*n_sig,ndof);
    eta = zeros(6*n_sig,1);
end


M_custom = zeros(ndof,ndof);
C_custom = zeros(ndof,ndof);
F_custom = zeros(ndof,1);

i_sig    = 1;
i_sig_nj = 1;

for i=1:N

    i_sig = i_sig+1; %joint

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        J_here    = J((i_sig-1)*6+1:i_sig*6,:);
        Jd_here   = Jd((i_sig-1)*6+1:i_sig*6,:);
        eta_here  = eta((i_sig-1)*6+1:i_sig*6);
        M_here    = M_added((i_sig_nj-1)*6+1:i_sig_nj*6,:);
        Fext_here = Fext((i_sig_nj-1)*6+1:i_sig_nj*6);
        M_custom  = M_custom+J_here'*M_here*J_here;
        C_custom  = C_custom+J_here'*(M_here*Jd_here+dinamico_coadj(eta_here)*M_here*J_here);
        F_custom  = F_custom+J_here'*Fext_here;
        i_sig_nj  = i_sig_nj+1;
        i_sig     = i_sig+1;
    end

    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        lpf        = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Ws         = Tr.VLinks(Tr.LinkIndex(i)).Ws{j};
        nGauss     = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};

        %scaling of quantities
        Lscale = lpf;
        lpf    = 1;

        i_sig      = i_sig+1;
        i_sig_nj   = i_sig_nj+1;
        for ii=2:nGauss-1
            J_here    = J((i_sig-1)*6+1:i_sig*6,:);
            Jd_here   = Jd((i_sig-1)*6+1:i_sig*6,:);
            eta_here  = eta((i_sig-1)*6+1:i_sig*6);
            M_here    = M_added((i_sig_nj-1)*6+1:i_sig_nj*6,:);
            Fext_here = Fext((i_sig_nj-1)*6+1:i_sig_nj*6);

            Fext_here(1:3)  = Fext_here(1:3)/Lscale;
            M_here(1:3,1:3) = M_here(1:3,1:3)/Lscale;
            M_here(4:6,4:6) = M_here(4:6,4:6)*Lscale;

            W_here   = Ws(ii);
            M_custom = M_custom+lpf*W_here*J_here'*M_here*J_here*Lscale^2; %scaling back
            F_custom = F_custom+lpf*W_here*J_here'*Fext_here*Lscale^2; %scaling back
            C_custom = C_custom+lpf*W_here*J_here'*(M_here*Jd_here+dinamico_coadj(eta_here)*M_here*J_here)*Lscale^2; %scaling back
            i_sig    = i_sig+1;
            i_sig_nj = i_sig_nj+1;
        end
        i_sig      = i_sig+1;
        i_sig_nj   = i_sig_nj+1;

    end
end

end

