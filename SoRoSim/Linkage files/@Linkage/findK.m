%Computes the generalized stiffness matrix for the linkage
%Last modified by Anup Teejo Mathew - 20/05/2021
function K = findK(S)

Ki = cell(1,S.ntot); 
f=1;
for i=1:S.N

    if S.VLinks(i).jointtype=='R'
        
        close all
        S.plotq0(f);
        prompt           = ['Enter the torsional stiffness of revolute joint of link ',num2str(i),' (Nm/rad):'];
        dlgtitle         = 'Torsional Stiffness';
        definput         = {'0'};
        opts.Interpreter = 'tex';
        opts.WindowStyle = 'normal';
        KR               = inputdlg(prompt,dlgtitle,[1 75],definput,opts);
        Ki{f}            = str2double(KR{1});
        
    elseif S.VLinks(i).jointtype=='P'
        
        close all
        S.plotq0(f);
        prompt           = ['Enter the linear stiffness of prismatic joint of link ',num2str(i),' (N/m):'];
        dlgtitle         = 'Linear Stiffness';
        definput         = {'0'};
        opts.Interpreter = 'tex';
        opts.WindowStyle = 'normal';
        KR               = inputdlg(prompt,dlgtitle,[1 75],definput,opts);
        Ki{f}            = str2double(KR{1});
        
    elseif S.VLinks(i).jointtype=='N'
        Ki{f} = [];
    else
        Ki{f} = zeros(S.Vtwists(f).dof);
    end
    
    %Add stiffness for more joints here
    
    f=f+1;
         
    for j=1:(S.VLinks(i).npie)-1 %for the rest of soft link divisions
        lpf     = S.VLinks(i).lp{j};
        Es     = S.VLinks(i).Es{j};
        B      = S.Vtwists(f).B;
        dof    = S.Vtwists(f).dof;
        Ktemp  = zeros(dof,dof);
        Ws     = S.VLinks(i).Ws{j};
        nGauss = S.VLinks(i).nGauss{j};


        %scaling of quantities
        Lscale        = lpf;
        lpf           = 1;
        
        for ii=1:nGauss
            Es_here = Es((ii-1)*6+1:ii*6,:);
            %scaling
            Es_here(1:3,:) = Es_here(1:3,:)/Lscale^3;
            Es_here(4:6,:) = Es_here(4:6,:)/Lscale;
            Ktemp   = Ktemp+lpf*Ws(ii)*B((ii-1)*6+1:ii*6,:)'*Es_here*B((ii-1)*6+1:ii*6,:);
        end
        Ki{f} = Ktemp*Lscale^2; %scaling back 
        
        f=f+1;
        
    end
    
end

%construct K matrix:
K = [];
for i=1:S.ntot
    K = blkdiag(K,Ki{i});
end
K=K./(S.q_scale*S.q_scale'); %actual K
end