%Function that allows the user to specify prismatic joint control and actuation
%specifications (24.05.2021)

function [n_Pact,i_Pact,i_Pactq,WrenchControlledP,BqP] = PrismaticJointActuation(Tr)
ndof              = Tr.ndof;
B1                = zeros(ndof,1);

n_Pact            = 0;
i_Pact            = [];
i_Pactq           = [];
BqP               = [];
WrenchControlledP = [];

dofi              = 1;

for i = 1:Tr.N %for each link
    
    VTwists_i = Tr.CVTwists{i};
    
    if Tr.VLinks(Tr.LinkIndex(i)).jointtype == 'P'
        
        close all
        Tr.plotq0(i);
        
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
                        WrenchControlledP = [WrenchControlledP;true];
                    case 'Displacement'
                        WrenchControlledP = [WrenchControlledP;false];
                end
                
                dofi = dofi+1;
        end
    else
        dofi = dofi+VTwists_i(1).dof;
    end
    
    for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        dofi = dofi+VTwists_i(j+1).dof;
    end
end

end
