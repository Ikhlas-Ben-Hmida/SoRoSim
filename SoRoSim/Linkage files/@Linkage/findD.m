%Computes the generalized damping matrix for the linkage
%Last modified by Anup Teejo Mathew - 23/05/2021
function D = findD(S)

Di = cell(1,S.ntot); 
f=1;
for i=1:S.N
    %for joint (rigid link)
    if S.VLinks(i).jointtype=='N'
        Di{f} = [];
    else
        Di{f} = zeros(S.Vtwists(f).dof);
    end
    
    f=f+1;  
    for j=1:(S.VLinks(i).npie)-1 %for the rest of soft pieces

        lpf    = S.VLinks(i).lp{j};
        Gs     = S.VLinks(i).Gs{j};
        B      = S.Vtwists(f).B;
        dof    = S.Vtwists(f).dof;
        Dtemp  = zeros(dof,dof);
        Ws     = S.VLinks(i).Ws{j};
        nGauss = S.VLinks(i).nGauss{j};

        %scaling of quantities
        Lscale        = lpf;
        lpf           = 1;

        for ii=1:nGauss
            Gs_here = Gs((ii-1)*6+1:ii*6,:);
            %scaling
            Gs_here(1:3,:) = Gs_here(1:3,:)/Lscale^3;
            Gs_here(4:6,:) = Gs_here(4:6,:)/Lscale;
            Dtemp = Dtemp+lpf*Ws(ii)*B((ii-1)*6+1:ii*6,:)'*Gs_here*B((ii-1)*6+1:ii*6,:);
        end

        Di{f} = Dtemp*Lscale^2; %scaling back 

        f=f+1;
    end
    
end

%construct K matrix:
D = [];
for i=1:S.ntot
    D = blkdiag(D,Di{i});
end
D=D./(S.q_scale*S.q_scale'); %actual K
end
%end
