%Function that allows the user to specify cylinderical joint control and actuation
%specifications (24.05.2021)

function [n_Cact,i_Cact,i_Cactq,WrenchControlledC] = CylindricalJointActuation(S)

n_Cact            = 0;
i_Cact            = [];
i_Cactq           = [];
WrenchControlledC = [];

f                 = 1;
dofi              = 1;

for i=1:S.N
    
    if S.VLinks(i).jointtype == 'C'
        
        close all
        S.plotq0(f);
        
        quest  = ['Is the cylindircal joint of link ',num2str(i),' actuated?'];
        answer = questdlg(quest,'Cylindrical Joint',...
            'Yes','No','Yes');
        
        switch answer
            case 'Yes'
                
                n_Cact  = n_Cact+2;
                i_Cact  = [i_Cact i];
                i_Cactq = [i_Cactq dofi dofi+1];
                
                quest = 'Is the cylindrical joint controlled by wrenches or joint coordinates (qs)?';
                answer2 = questdlg(quest,'Control',...
                    'Wrench','q','Wrench');
                
                switch answer2
                    case 'Wrench'
                        WrenchControlledC = [WrenchControlledC 1 1];
                    case 'q'
                        WrenchControlledC = [WrenchControlledC 0 0];
                end
                
                dofi    = dofi+2;
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
