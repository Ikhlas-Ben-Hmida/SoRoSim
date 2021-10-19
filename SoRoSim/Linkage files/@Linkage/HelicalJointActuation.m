%Function that allows the user to specify helical joint control and actuation
%specifications (24.05.2021)

function [n_Hact,i_Hact,i_Hactq,WrenchControlledH,BqH] = HelicalJointActuation(S)
ndof              = S.ndof;
B1                = zeros(ndof,1);

n_Hact            = 0;
i_Hact            = [];
i_Hactq           = [];
BqH               = [];
WrenchControlledH = [];

f                 = 1;
dofi              = 1;

for i = 1:S.N %for each link
    
    if S.VLinks(i).jointtype == 'H'
        
        close all
        S.plotq0(f);
        
        quest  = ['Is the Helical joint of link ',num2str(i),' actuated?'];
        answer = questdlg(quest,'Helical Joint',...
            'Yes','No','Yes');
        
        switch answer
            case 'Yes'
                
                n_Hact   = n_Hact+1;
                i_Hact   = [i_Hact;i];
                i_Hactq  = [i_Hactq;dofi];
                B1(dofi) = 1;
                BqH      = [BqH,B1];
                B1       = zeros(ndof,1);
                
                quest    = 'Is the helical joint controlled by wrench or joint coordinate (q)?';
                answer2  = questdlg(quest,'Control',...
                    'Wrench','q','Wrench');
                
                switch answer2
                    case 'Wrench'
                        WrenchControlledH = [WrenchControlledH;1];
                    case 'q'
                        WrenchControlledH = [WrenchControlledH;0];
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
