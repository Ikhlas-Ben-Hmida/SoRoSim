%Finds the scaling factor of joint coordinates
%Last modified by Anup Teejo Mathew 02.03.2022
function q_scale = find_q_scale(Tr)

ndof = 0;
for i=1:Tr.N
    ndof = ndof+Tr.CVTwists{i}(1).dof;
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        ndof = ndof+Tr.CVTwists{i}(j+1).dof;
    end
end
q_scale  = ones(ndof,1);

f=1;

for i=1:Tr.N %for each link
                   
    VTwists = Tr.CVTwists{i};
    f       = f+VTwists(1).dof;
    
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %for each of the soft link divisions
        
        Lscale = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        for ii=1:6
            if ii<4
                if VTwists(j+1).Bdof(ii)
                    for kk=1:VTwists(j+1).Bodr(ii)+1
                        q_scale(f) = 1/Lscale^kk;
                        f=f+1;
                    end
                end
            else
                if VTwists(j+1).Bdof(ii)
                    for kk=1:VTwists(j+1).Bodr(ii)+1
                        q_scale(f) = 1/Lscale^(kk-1);
                        f=f+1;
                    end
                end
            end
        end
        
    end

end

end

