%Computes the generalized damping matrix for the linkage
%Last modified by Anup Teejo Mathew 02.03.2022
function D = findD(Tr)

CD = cell(Tr.N,1); 

for i=1:Tr.N
    
    CD{i}   = cell(Tr.VLinks(Tr.LinkIndex(i)).npie,1); 
    VTwists = Tr.CVTwists{i};
    %for joint (rigid link)
    if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='N'
        CD{i}{1} = [];
    else
        CD{i}{1} = zeros(VTwists(1).dof);
    end
 
    for j=1:(Tr.VLinks(Tr.LinkIndex(i)).npie)-1 %for the rest of soft pieces

        ld     = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Gs     = Tr.VLinks(Tr.LinkIndex(i)).Gs{j};
        B      = VTwists(j+1).B;
        dof    = VTwists(j+1).dof;
        Dtemp  = zeros(dof,dof);
        Ws     = Tr.VLinks(Tr.LinkIndex(i)).Ws{j};
        nGauss = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};

        %scaling of quantities
        Lscale        = ld;
        ld           = 1;

        for ii=1:nGauss
            Gs_here = Gs((ii-1)*6+1:ii*6,:);
            %scaling
            Gs_here(1:3,:) = Gs_here(1:3,:)/Lscale^3;
            Gs_here(4:6,:) = Gs_here(4:6,:)/Lscale;
            
            Dtemp = Dtemp+ld*Ws(ii)*B((ii-1)*6+1:ii*6,:)'*Gs_here*B((ii-1)*6+1:ii*6,:);
        end
        CD{i}{j+1} = Dtemp*Lscale^2; %scaling back 

    end
    
end

%construct D matrix:
D = [];
for i=1:Tr.N
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie
        D = blkdiag(D,CD{i}{j});
    end
end

q_scale = find_q_scale(Tr);
D       = D./(q_scale*q_scale'); %actual D
end
