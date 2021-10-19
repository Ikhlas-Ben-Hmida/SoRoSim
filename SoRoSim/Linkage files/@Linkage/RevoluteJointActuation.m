%Function that allows the user to specify revolute joint control and actuation
%specifications (24.05.2021)

function [n_Ract,i_Ract,i_Ractq,WrenchControlledR,BqR] = RevoluteJointActuation(S)
ndof              = S.ndof;
B1                = zeros(ndof,1);

n_Ract            = 0;
i_Ract            = [];
i_Ractq           = [];
BqR               = [];
WrenchControlledR = [];

f                 = 1;
dofi              = 1;

for i=1:S.N %for each link
    
    if S.VLinks(i).jointtype=='R'
        
        close all
        S.plotq0(f);
        quest  = ['Is the revolute joint of link ',num2str(i),' actuated?'];
        answer = questdlg(quest,'Revolute Joint',...
            'Yes','No','Yes');
        
        switch answer
            case 'Yes'
                
                n_Ract   = n_Ract+1;
                i_Ract   = [i_Ract;i];
                i_Ractq  = [i_Ractq;dofi];
                B1(dofi) = 1;
                BqR      = [BqR,B1];
                B1       = zeros(ndof,1);
                
                quest    = 'Is the revolute joint controlled by torque or angle?';
                answer2  = questdlg(quest,'Control',...
                    'Torque','Angle','Torque');
                
                switch answer2
                    case 'Torque'
                        WrenchControlledR = [WrenchControlledR;1];
                    case 'Angle'
                        WrenchControlledR = [WrenchControlledR;0];
                end
                
                dofi = dofi+1;
        end
    else
        dofi = dofi+S.Vtwists(f).dof;
    end
    f = f+1;
    
    for j = 1:S.VLinks(i).npie-1
        dofi = dofi+S.Vtwists(f).dof;
        f    = f+1;
    end
end

end
