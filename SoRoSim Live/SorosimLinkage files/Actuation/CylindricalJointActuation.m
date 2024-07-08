%Function that allows the user to specify cylinderical joint control and actuation
%specifications (24.05.2021)

function [n_Cact,i_Cact,i_Cactq,WrenchControlledC] = CylindricalJointActuation(Tr,Update)
if nargin==1
    Update=false;
end
n_Cact            = 0;
i_Cact            = [];
i_Cactq           = [];
WrenchControlledC = [];

dofi              = 1;

if ~Update
    for i=1:Tr.N

        VTwists_i = Tr.CVTwists{i};

        if Tr.VLinks(Tr.LinkIndex(i)).jointtype == 'C'

            close all
            Tr.plotq0(i);

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
                            WrenchControlledC = [WrenchControlledC true true];
                        case 'q'
                            WrenchControlledC = [WrenchControlledC false false];
                    end


            end
        end
        dofi = dofi+VTwists_i(1).dof;
        for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
            dofi = dofi+VTwists_i(j+1).dof;
        end
    end
else
    for i=1:Tr.N %for each link

        VTwists_i = Tr.CVTwists{i};

        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='C'&&any(Tr.i_jact==i)
            
            i_Cactq  = [i_Cactq dofi dofi+1];

        end

        dofi = dofi+VTwists_i(1).dof;
        for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
            dofi = dofi+VTwists_i(j+1).dof;
        end
    end
end

end
