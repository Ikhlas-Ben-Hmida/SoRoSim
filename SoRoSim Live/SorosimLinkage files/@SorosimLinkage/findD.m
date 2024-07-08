%Computes the generalized damping matrix for the linkage
%Last modified by Anup Teejo Mathew 14.04.2023
function D = findD(Tr,varargin)

dof_start = 1;

D = zeros(Tr.ndof,Tr.ndof); %NxN matrix 
for i=1:Tr.N

    VTwists = Tr.CVTwists{i};
    dof_here = VTwists(1).dof;
    %for joint (rigid link)
    if Tr.VLinks(Tr.LinkIndex(i)).Dj==0
        D(dof_start:dof_start+dof_here-1,dof_start:dof_start+dof_here-1) = zeros(VTwists(1).dof);
    else
        if isequal(size(Tr.VLinks(Tr.LinkIndex(i)).Dj),[VTwists(1).dof,VTwists(1).dof])
            D(dof_start:dof_start+dof_here-1,dof_start:dof_start+dof_here-1) = Tr.VLinks(Tr.LinkIndex(i)).Dj;
        else
            uiwait(msgbox('Incorrect joint damping matrix dimensions','Error','error'));
            return
        end
    end
    dof_start = dof_start+dof_here;
    for j=1:(Tr.VLinks(Tr.LinkIndex(i)).npie)-1 %for the rest of soft pieces

        ld       = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Gs       = VTwists(j+1).Gs;
        dof_here = VTwists(j+1).dof;
        Ws       = VTwists(j+1).Ws;
        nip       = VTwists(j+1).nip;

        Dtemp  = zeros(dof_here,dof_here);

        %scaling of quantities
        Lscale = ld;
        dBqdq  = VTwists(j+1).B;

        for ii=1:nip
            if Ws(ii)>0
                Gs_here = Gs((ii-1)*6+1:ii*6,:);
                %scaling
                Gs_here(1:3,:) = Gs_here(1:3,:)/Lscale^3;
                Gs_here(4:6,:) = Gs_here(4:6,:)/Lscale;

                Dtemp = Dtemp+Ws(ii)*dBqdq((ii-1)*6+1:ii*6,:)'*Gs_here*dBqdq((ii-1)*6+1:ii*6,:);
            end
        end

        D(dof_start:dof_start+dof_here-1,dof_start:dof_start+dof_here-1) = Dtemp*Lscale^2; %scaling back 
        dof_start  = dof_start+dof_here;
    end
end

end
