%Computes the generalized stiffness matrix for the linkage
%Last modified by Anup Teejo Mathew 02.03.2022
function K = findK(Tr)

CK = cell(Tr.N,1);

for i=1:Tr.N
    
    CK{i}   = cell(Tr.VLinks(Tr.LinkIndex(i)).npie,1);
    VTwists = Tr.CVTwists{i};

    if Tr.VLinks(Tr.LinkIndex(i)).Kj==0
        CK{i}{1} = zeros(VTwists(1).dof);
    else
        if isequal(size(Tr.VLinks(Tr.LinkIndex(i)).Kj),[VTwists(1).dof,VTwists(1).dof])
            CK{i}{1} = Tr.VLinks(Tr.LinkIndex(i)).Kj;
        else
            uiwait(msgbox('Incorrect stiffness matrix dimensions','Error','error'));
            return
        end
    end
    
    %Add stiffness for more joints here
         
    for j=1:(Tr.VLinks(Tr.LinkIndex(i)).npie)-1 %for the rest of soft link divisions
        
        ld     = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        Es     = Tr.VLinks(Tr.LinkIndex(i)).Es{j};
        B      = VTwists(j+1).B;
        dof    = VTwists(j+1).dof;
        Ktemp  = zeros(dof,dof);
        Ws     = Tr.VLinks(Tr.LinkIndex(i)).Ws{j};
        nGauss = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};

        %scaling of quantities
        Lscale        = ld;
        ld            = 1;
        
        for ii=1:nGauss
            Es_here = Es((ii-1)*6+1:ii*6,:);
            %scaling
            Es_here(1:3,:) = Es_here(1:3,:)/Lscale^3;
            Es_here(4:6,:) = Es_here(4:6,:)/Lscale;
            
            Ktemp = Ktemp+ld*Ws(ii)*B((ii-1)*6+1:ii*6,:)'*Es_here*B((ii-1)*6+1:ii*6,:);
        end
        CK{i}{j+1} = Ktemp*Lscale^2; %scaling back 
        
    end
    
end

%construct K matrix:
K = [];
for i=1:Tr.N
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie
        K = blkdiag(K,CK{i}{j});
    end
end

q_scale = find_q_scale(Tr);
K       = K./(q_scale*q_scale'); %actual K
end