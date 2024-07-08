%Function that allows the user to specify revolute joint control and actuation
%specifications (24.05.2021)

function [n_Ract,i_Ract,i_Ractq,WrenchControlledR,BqR] = RevoluteJointActuation(Tr,Update)
if nargin==1
    Update=false;
end
ndof              = Tr.ndof;
B1                = zeros(ndof,1);

n_Ract            = 0;
i_Ract            = [];
i_Ractq           = [];
BqR               = [];
WrenchControlledR = [];

dofi              = 1;
if ~Update
    for i=1:Tr.N %for each link

        VTwists_i = Tr.CVTwists{i};

        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='R'

            close all
            Tr.plotq0(i);
            quest  = ['Is the revolute joint of link ',num2str(i),' actuated?'];
            answer = questdlg(quest,'Revolute Joint',...
                'Yes','No','Yes');

            switch answer
                case 'Yes'

                    n_Ract   = n_Ract+1;
                    i_Ract   = [i_Ract i];
                    i_Ractq  = [i_Ractq dofi];
                    B1(dofi) = 1;
                    BqR      = [BqR,B1];
                    B1       = zeros(ndof,1);

                    quest    = 'Is the revolute joint controlled by torque or angle?';
                    answer2  = questdlg(quest,'Control',...
                        'Torque','Angle','Torque');

                    switch answer2
                        case 'Torque'
                            WrenchControlledR = [WrenchControlledR true];
                        case 'Angle'
                            WrenchControlledR = [WrenchControlledR false];
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

        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='R'&&any(Tr.i_jact==i)
            
            i_Ractq  = [i_Ractq dofi];
            B1(dofi) = 1;
            BqR      = [BqR,B1];
            B1       = zeros(ndof,1);

        end

        dofi = dofi+VTwists_i(1).dof;
        for j = 1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
            dofi = dofi+VTwists_i(j+1).dof;
        end
    end
end

end
