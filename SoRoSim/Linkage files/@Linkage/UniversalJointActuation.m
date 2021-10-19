%Function that allows the user to specify universal joint control and actuation
%specifications (24.05.2021)

function [n_Uact,i_Uact,i_Uactq,WrenchControlledU] = UniversalJointActuation(S)

n_Uact            = 0;
i_Uact            = [];
i_Uactq           = [];
WrenchControlledU = [];

f                 = 1;
dofi              = 1;

for i = 1:S.N %for each link
    
    if S.VLinks(i).jointtype == 'U'
        
        close all
        S.plotq0(f);
        
        quest  = ['Is the universal joint of link ',num2str(i),' actuated?'];
        answer = questdlg(quest,'Universal Joint',...
            'Yes','No','Yes');
        
        switch answer
            case 'Yes'
                
                n_Uact  = n_Uact+2;
                i_Uact  = [i_Uact i];
                i_Uactq = [i_Uactq dofi dofi+1];
                
                quest   = 'Is the universal joint controlled by wrenches or joint coordinates (qs)?';
                answer2 = questdlg(quest,'Control',...
                    'Wrench','q','Wrench');
                
                switch answer2
                    case 'Wrench'
                        WrenchControlledU = [WrenchControlledU 1 1];
                    case 'q'
                        WrenchControlledU = [WrenchControlledU 0 0];
                end
                
                dofi = dofi+2;
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
