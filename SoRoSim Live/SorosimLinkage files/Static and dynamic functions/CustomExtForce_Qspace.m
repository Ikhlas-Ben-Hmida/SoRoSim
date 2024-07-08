%Function to convert user defined external force (Fext) on each significant
%points to external force in Q space
%Last modified by Anup Teejo Mathew - 24/05/2021
function F_custom=CustomExtForce_Qspace(Tr,J,Fext,FextP)

N        = Tr.N;
ndof     = Tr.ndof;

F_custom = zeros(ndof,1);

i_sig    = 1; %index of all significant points
i_sig_nj = 1; %significant points excluding joints

for i=1:N

    i_sig = i_sig+1; %joint

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        J_here    = J((i_sig-1)*6+1:i_sig*6,:);
        Fext_here = Fext((i_sig_nj-1)*6+1:i_sig_nj*6)+FextP((i_sig_nj-1)*6+1:i_sig_nj*6);
        F_custom  = F_custom+J_here'*Fext_here;
        i_sig_nj  = i_sig_nj+1;
        i_sig     = i_sig+1;
    end

    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        
        lpf        = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Ws         = Tr.CVTwists{i}(j+1).Ws;
        nip        = Tr.CVTwists{i}(j+1).nip;

        %scaling of quantities
        Lscale = lpf;
        for ii=1:nip
            W_here        = Ws(ii);
            J_here         = J((i_sig-1)*6+1:i_sig*6,:);
            Fext_here      = Fext((i_sig_nj-1)*6+1:i_sig_nj*6);
            Fext_here(1:3) = Fext_here(1:3)/Lscale;
            
            FextP_here      = FextP((i_sig_nj-1)*6+1:i_sig_nj*6);
            FextP_here(1:3) = FextP_here(1:3)/Lscale^2;
            FextP_here(4:6) = FextP_here(4:6)/Lscale;
            
            F_custom = F_custom+W_here*J_here'*Fext_here*Lscale^2 + J_here'*FextP_here*Lscale^2; %scaling back
            i_sig    = i_sig+1;
            i_sig_nj = i_sig_nj+1;
        end

    end
end

end

