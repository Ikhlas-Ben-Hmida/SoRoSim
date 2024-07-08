%Function that calculates the screw velocity (eta) at every significant point
%in SI units (last modified 18/04/2023, Anup)
function xi = ScrewStrain(Tr,q,i_here,j_here) %i_here is link index, j_here is division (0 for joint)

if isrow(q)
    q=q';
end

N         = Tr.N;

full = false;
if nargin==4
    if j_here==0
        nsig = 1;
    else
        nsig = Tr.CVTwists{i_here}(j_here+1).nip; %j_here>1 is allowed only for soft links
    end
elseif nargin==3
    if Tr.VLinks(Tr.LinkIndex(i_here)).linktype=='s'
        nsig = 1;
        for j=1:Tr.VLinks(Tr.LinkIndex(i_here)).npie-1
            nsig = nsig+Tr.CVTwists{i_here}(j+1).nip;
        end
    else
        nsig = 2; %joint and CM
    end
else
    full  = true;
    nsig  = Tr.nsig;
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
        xi_here = zeros(6,1);
    else
        xi_here                   = B_here*q_here+xi_star;
    end
    if full||(i==i_here&&nargin==3)||(i==i_here&&j_here==0)
        xi((i_sig-1)*6+1:i_sig*6) = xi_here;
        i_sig                     = i_sig+1;
    end
    
    if Tr.VLinks(Tr.LinkIndex(i)).linktype == 'r'
        
        if full||(i==i_here&&nargin==3)
            xi((i_sig-1)*6+1:i_sig*6) = zeros(6,1); % doesnt have a meaning
            i_sig                     = i_sig+1;
        end

    end
    
    dof_start = dof_start+dof_here;
    
    for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        
        dof_here = Tr.CVTwists{i}(j+1).dof;
        q_here   = q(dof_start:dof_start+dof_here-1);
        xi_star  = Tr.CVTwists{i}(j+1).xi_star;
        Lscale   = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        nip      = Tr.CVTwists{i}(j+1).nip;
        B        = Tr.CVTwists{i}(j+1).B;

        
        for ii = 1:nip
            
             xi_here = xi_star(6*(ii-1)+1:6*ii,1);
             xi_here(1:3) = xi_here(1:3)*Lscale; %scaling using the formula: Lscale m = 1 unit

            B_here  = B(6*(ii-1)+1:6*ii,:);%note this step
            if dof_here>0
                xi_here = B_here*q_here+xi_here; 
            end
            
            xi_here(1:3) = xi_here(1:3)/Lscale; %to SI units
            
            if full||(i==i_here&&nargin==3)||(i==i_here&&j==j_here)
                xi((i_sig-1)*6+1:i_sig*6) = xi_here;
                i_sig                     = i_sig+1;
            end
            
        end
        
        dof_start  = dof_start+dof_here;
    end
    
end
end

