%Finds the local initial transformation matrix for a link towards center of
%gravity for rigid link, towards center of area for soft link
%Last modified by Anup Teejo Mathew - 23/05/2021
function g0 = findg0(S)
%=====================================================================
% initial transformation (g0)

g0p = eye(4); %initializing g0
g0  = cell(1,S.ntot);

f=1; %index of piece
for i=1:S.N
    
        if S.VLinks(i).linktype=='r' %if rigid link
            g0p(:,4) = [S.VLinks(i).cx;S.VLinks(i).cy;S.VLinks(i).cz;1];
            g0{f}    = g0p;
            g0p      = eye(4);
            
            f=f+1;
        elseif S.VLinks(i).linktype=='s' %if soft link
            for j=1:S.VLinks(i).npie
                if j==1
                    g0{f} = eye(4);%for joint
                    
                    f=f+1;
                else
                    g0p(:,4) = [0;S.VLinks(i).cy{j-1};S.VLinks(i).cz{j-1};1];
                    g0{f}    = g0p;
                    g0p      = eye(4);
                    
                    f=f+1;
                end
            end
        end
        
end

