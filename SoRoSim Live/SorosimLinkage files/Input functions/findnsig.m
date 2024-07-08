%Finds the number of points at which computations are performed (significant points)
%Last modified by Anup Teejo Mathew 02.03.2022
function nsig = findnsig(Tr)

nsig = Tr.N; %all joints
for i=1:Tr.N
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        nsig = nsig+Tr.CVTwists{i}(j+1).nip;
    end
    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        nsig = nsig+1;
    end
end

end