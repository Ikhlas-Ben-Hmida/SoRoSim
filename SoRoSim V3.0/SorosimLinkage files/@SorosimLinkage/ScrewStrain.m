%Function that calculates the screw velocity (eta) at every significant point
%(24.05.2021)

function xi = ScrewStrain(Tr,q)

if isrow(q)
    q=q';
end

N         = Tr.N;
nsig      = Tr.nsig;

%% Mass, Corriolis, Gravity
for i=1:N
    if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='U' %two strain values must be saved for a universal joint.
        nsig = nsig+1; 
    end
end

xi        = zeros(6*nsig,1);
dof_start = 1; %starting dof of current piece
i_sig     = 1;

for i = 1:N
    %Joint
    dof_here = Tr.CVTwists{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    B_here   = Tr.CVTwists{i}(1).B;
    xi_star  = Tr.CVTwists{i}(1).xi_star;
    
    if dof_here == 0 %fixed joint (N)

    else
        if Tr.VLinks(Tr.LinkIndex(i)).jointtype == 'U' %special case for universal joint. Considered as 2 revolute joints
            % first revolute joint
            xi_here                   = B_here(:,1)*q_here(1)+xi_star;
            xi((i_sig-1)*6+1:i_sig*6) = xi_here;
            i_sig                     = i_sig+1;
                     
            % second revolute joint
            xi_here                   = B_here(:,2)*q_here(2)+xi_star;
            xi((i_sig-1)*6+1:i_sig*6) = xi_here;
            i_sig                     = i_sig+1;
            
        else
            xi_here                   = B_here*q_here+xi_star;
            xi((i_sig-1)*6+1:i_sig*6) = xi_here;
            i_sig                     = i_sig+1;
        end
    end
    
    if Tr.VLinks(Tr.LinkIndex(i)).linktype == 'r'
        
        xi_sig((i_sig-1)*6+1:i_sig*6)    = zeros(6,1);
        i_sig                            = i_sig+1;

    end
    
    dof_start = dof_start+dof_here;
    
    for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        
        dof_here = Tr.CVTwists{i}(j+1).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);
        Bdof     = Tr.CVTwists{i}(j+1).Bdof;
        Bodr     = Tr.CVTwists{i}(j+1).Bodr;
        xi_star  = Tr.CVTwists{i}(j+1).xi_star;
        nGauss   = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};
        lpf      = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        
        q_scale_here  = Tr.q_scale(dof_start:dof_start+dof_here-1);
        doftheta_here = Bdof(1:3)'*(Bodr(1:3)+[1 1 1]');
        q_scale_here(1:doftheta_here) = q_scale_here(1:doftheta_here)*lpf;
        B_scale = repmat(q_scale_here',6*nGauss,1);
        B       = Tr.CVTwists{i}(j+1).B./B_scale; %actual 
        
        for ii = 1:nGauss
            
            B_here      = B(6*(ii-1)+1:6*ii,:);
            xi_starhere = xi_star(6*(ii-1)+1:6*ii,1);

            xi_here                   = B_here*q_here+xi_starhere;
            xi((i_sig-1)*6+1:i_sig*6) = xi_here;
            i_sig                     = i_sig+1;
            
        end
        %updating g, Jacobian, Jacobian_dot and eta at X=L
        dof_start  = dof_start+dof_here;
    end
    
end
end

