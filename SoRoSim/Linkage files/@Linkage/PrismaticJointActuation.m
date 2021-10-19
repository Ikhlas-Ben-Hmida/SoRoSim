%Function that allows the user to specify prismatic joint control and actuation
%specifications (24.05.2021)

function [n_Pact,i_Pact,i_Pactq,WrenchControlledP,BqP] = PrismaticJointActuation(S)
ndof              = S.ndof;
B1                = zeros(ndof,1);

n_Pact            = 0;
i_Pact            = [];
i_Pactq           = [];
BqP               = [];
WrenchControlledP = [];

f                 = 1;
dofi              = 1;

for i = 1:S.N %for each link
    
    if S.VLinks(i).jointtype == 'P'
        
        close all
        S.plotq0(f);
        
        quest  = ['Is the prismatic joint of link ',num2str(i),' actuated?'];
        answer = questdlg(quest,'Prismatic Joint',...
            'Yes','No','Yes');
        
        switch answer
            case 'Yes'
                
                n_Pact   = n_Pact+1;
                i_Pact   = [i_Pact;i];
                i_Pactq  = [i_Pactq;dofi];
                B1(dofi) = 1;
                BqP      = [BqP,B1];
                B1       = zeros(ndof,1);
                
                quest    = 'Is the prismatic joint controlled by force or displacement?';
                answer2  = questdlg(quest,'Control',...
                    'Force','Displacement','Force');
                
                switch answer2
                    case 'Force'
                        WrenchControlledP = [WrenchControlledP;1];
                    case 'Displacement'
                        WrenchControlledP = [WrenchControlledP;0];
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
