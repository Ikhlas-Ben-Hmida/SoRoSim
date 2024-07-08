%Computes the generalized stiffness matrix for the linkage
%Last modified by Anup Teejo Mathew 14.04.2023
function K = findK(Tr,varargin) %NxN matrix for independent basis, NX1 for dependent basis if any

dof_start = 1;

K = zeros(Tr.ndof,Tr.ndof); %NxN matrix 
for i=1:Tr.N

    VTwists = Tr.CVTwists{i};
    dof_here = VTwists(1).dof;
    if Tr.VLinks(Tr.LinkIndex(i)).Kj==0
        K(dof_start:dof_start+dof_here-1,dof_start:dof_start+dof_here-1) = zeros(VTwists(1).dof);
    else
        if isequal(size(Tr.VLinks(Tr.LinkIndex(i)).Kj),[VTwists(1).dof,VTwists(1).dof])
            K(dof_start:dof_start+dof_here-1,dof_start:dof_start+dof_here-1) = Tr.VLinks(Tr.LinkIndex(i)).Kj;
        else
            uiwait(msgbox('Incorrect joint stiffness matrix dimensions','Error','error'));
            return
        end
    end

    dof_start = dof_start+dof_here;  
    for j=1:(Tr.VLinks(Tr.LinkIndex(i)).npie)-1 %for the rest of soft link divisions

        ld       = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Es       = VTwists(j+1).Es;
        dof_here = VTwists(j+1).dof;
        Ws       = VTwists(j+1).Ws;
        nip       = VTwists(j+1).nip;

        Ktemp  = zeros(dof_here,dof_here);

        %scaling of quantities
        Lscale = ld;
        dBqdq  = VTwists(j+1).B;

        for ii=1:nip
            if Ws(ii)>0
                Es_here = Es((ii-1)*6+1:ii*6,:);
                %scaling
                Es_here(1:3,:) = Es_here(1:3,:)/Lscale^3;
                Es_here(4:6,:) = Es_here(4:6,:)/Lscale;

                Ktemp = Ktemp+Ws(ii)*dBqdq((ii-1)*6+1:ii*6,:)'*Es_here*dBqdq((ii-1)*6+1:ii*6,:);
            end
        end
        K(dof_start:dof_start+dof_here-1,dof_start:dof_start+dof_here-1) = Ktemp*Lscale^2;
        dof_start  = dof_start+dof_here;
    end
end

end