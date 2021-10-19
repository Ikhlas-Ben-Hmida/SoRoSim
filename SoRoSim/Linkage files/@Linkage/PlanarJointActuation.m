%Function that allows the user to specify planar joint control and actuation
%specifications (24.05.2021)

function [n_Aact,i_Aact,i_Aactq,WrenchControlledA] = PlanarJointActuation(S)

n_Aact            = 0;
i_Aact            = [];
i_Aactq           = [];
WrenchControlledA = [];

f                 = 1;
dofi              = 1;

for i=1:S.N %for each link
    
    if S.VLinks(i).jointtype == 'A'
        
        close all
        S.plotq0(f);
        
        quest  = ['Is the planar joint of link ',num2str(i),' actuated?'];
        answer = questdlg(quest,'Planar Joint',...
            'Yes','No','Yes');
        
        switch answer
            case 'Yes'
                
                n_Aact  = n_Aact+3;
                i_Aact  = [i_Aact i];
                i_Aactq = [i_Aactq dofi dofi+1 dofi+2];
                
                quest   = 'Is the planar joint controlled by wrenches or joint coordinates (qs)?';
                answer2 = questdlg(quest,'Control',...
                    'Wrench','q','Wrench');
                
                switch answer2
                    case 'Wrench'
                        WrenchControlledA = [WrenchControlledA 1 1 1];
                    case 'q'
                        WrenchControlledA = [WrenchControlledA 0 0 0];
                end
                
                dofi = dofi+3;
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
