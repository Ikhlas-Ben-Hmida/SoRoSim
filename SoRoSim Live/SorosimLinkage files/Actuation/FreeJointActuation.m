%Function that allows the user to specify free joint control and actuation
%specifications (24.05.2021)

function [n_Fact,i_Fact,i_Factq,WrenchControlledF] = FreeJointActuation(Tr,Update)
if nargin==1
    Update=false;
end
n_Fact            = 0;
i_Fact            = [];
i_Factq           = [];
WrenchControlledF = [];
dofi              = 1;
if ~Update
    for i=1:Tr.N %for each link

        VTwists_i = Tr.CVTwists{i};

        if Tr.VLinks(Tr.LinkIndex(i)).jointtype == 'F'

            close all
            Tr.plotq0(i);

            quest  = ['Is the Free joint of link ',num2str(i),' actuated?'];
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
                            WrenchControlledF = [WrenchControlledF true true true true true true];
                        case 'q'
                            WrenchControlledF = [WrenchControlledF false false false false false false];
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

        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='F'&&any(Tr.i_jact==i)
            
            i_Factq = [i_Factq dofi dofi+1 dofi+2 dofi+3 dofi+4 dofi+5];

        end

        dofi = dofi+VTwists_i(1).dof;
        for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
            dofi = dofi+VTwists_i(j+1).dof;
        end
    end
end
end
