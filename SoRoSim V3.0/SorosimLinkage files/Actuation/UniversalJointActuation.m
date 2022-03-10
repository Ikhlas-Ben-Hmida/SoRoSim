%Function that allows the user to specify universal joint control and actuation
%specifications (24.05.2021)

function [n_Uact,i_Uact,i_Uactq,WrenchControlledU] = UniversalJointActuation(Tr)

n_Uact            = 0;
i_Uact            = [];
i_Uactq           = [];
WrenchControlledU = [];

dofi              = 1;

for i = 1:Tr.N %for each link
    
    VTwists_i = Tr.CVTwists{i};
    
    if Tr.VLinks(Tr.LinkIndex(i)).jointtype == 'U'
        
        close all
        Tr.plotq0(i);
        
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
                        WrenchControlledU = [WrenchControlledU true true];
                    case 'q'
                        WrenchControlledU = [WrenchControlledU false false];
                end
        end
        
    end
    
    dofi = dofi+VTwists_i(1).dof;
    for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        dofi = dofi+VTwists_i(j+1).dof;
    end
end

end
