%Function that allows the user to specify spherical joint control and actuation
%specifications (24.05.2021)

function [n_Sact,i_Sact,i_Sactq,WrenchControlledS] = SphericalJointActuation(S)

n_Sact            = 0;
i_Sact            = [];
i_Sactq           = [];
WrenchControlledS = [];

f                 = 1;
dofi              = 1;

for i=1:S.N %for each link
    
    if S.VLinks(i).jointtype == 'S'
        
        close all
        S.plotq0(f);
        
        quest  = ['Is the spherical joint of link ',num2str(i),' actuated?'];
        answer = questdlg(quest,'Spherical Joint',...
            'Yes','No','Yes');
        
        switch answer
            case 'Yes'
                
                n_Sact  = n_Sact+3;
                i_Sact  = [i_Sact i];
                i_Sactq = [i_Sactq dofi dofi+1 dofi+2];
                
                quest   = 'Is the spherical joint controlled by wrenches or joint coordinates (qs)?';
                answer2 = questdlg(quest,'Control',...
                    'Wrench','q','Wrench');
                
                switch answer2
                    case 'Wrench'
                        WrenchControlledS = [WrenchControlledS 1 1 1];
                    case 'q'
                        WrenchControlledS = [WrenchControlledS 0 0 0];
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
