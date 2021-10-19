%Function that allows the user to specify free joint control and actuation
%specifications (24.05.2021)

function [n_Fact,i_Fact,i_Factq,WrenchControlledF] = FreeJointActuation(S)

n_Fact            = 0;
i_Fact            = [];
i_Factq           = [];
WrenchControlledF = [];

f                 = 1;
dofi              = 1;

for i=1:S.N %for each link
    
    if S.VLinks(i).jointtype == 'F'
        
        close all
        S.plotq0(f);
        
        quest  = ['Is the spherical joint of link ',num2str(i),' actuated?'];
        answer = questdlg(quest,'Spherical Joint',...
            'Yes','No','Yes');
        
        switch answer
            case 'Yes'
                
                n_Fact  = n_Fact+6;
                i_Fact  = [i_Fact i];
                i_Factq = [i_Factq dofi dofi+1 dofi+2 dofi+3 dofi+4 dofi+5];
                
                quest   = 'Is the spherical joint controlled by wrenches or joint coordinates (qs)?';
                answer2 = questdlg(quest,'Control',...
                    'Wrench','q','Wrench');
                
                switch answer2
                    case 'Wrench'
                        WrenchControlledF = [WrenchControlledF 1 1 1 1 1 1];
                    case 'q'
                        WrenchControlledF = [WrenchControlledF 0 0 0 0 0 0];
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
