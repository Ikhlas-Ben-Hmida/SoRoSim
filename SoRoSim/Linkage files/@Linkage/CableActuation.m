%Function that allows the user to set actuation cable parameters
%(24.05.2021)
function [nact,dc,dcp,Sdiv,Ediv,Inside] = CableActuation(S)

prompt           = {'Number of actuators:'};
dlgtitle         = 'Actuation Parameters';
definput         = {'1'};
opts.Interpreter = 'tex';
ans_act          = inputdlg(prompt,dlgtitle,[1 50],definput,opts);
nact             = str2num(ans_act{1});
N                = S.N;

dc               = cell(nact,N);
dcp              = cell(nact,N);
Sdiv             = cell(nact,N);
Ediv             = cell(nact,N);
Inside           = cell(1,nact);

%Set cable parameters
for ii = 1:nact
    
    if any(S.CS=='R')
        quest  = ['Is the actutator (cable) ',num2str(ii),' fully inside the linkage?'];
        Answer = questdlg(quest,'Inside?','Yes','No','Yes');
    else
        Answer = 'Yes';
    end

    switch Answer
        case 'Yes'
            Inside{ii} = 1;
            for i=1:N
                if S.VLinks(i).linktype=='s'
                    cableactuationGUI(S,ii,i)
                    load('cableactuation.mat','dci','dcpi','Sdivi','Edivi')
                    dc{ii,i}     = dci;
                    dcp{ii,i}    = dcpi;
                    Ediv{ii,i}   = Edivi;
                    Sdiv{ii,i}   = Sdivi;
                end
            end
        case 'No'
            Inside{ii} = 0;
            for i=1:N
                if S.VLinks(i).linktype=='s'
                    cableactuationGUI2(S,ii,i)
                    load('cableactuation.mat','dci','Sdivi','Edivi')
                    dc{ii,i}     = dci;
                    Ediv{ii,i}   = Edivi;
                    Sdiv{ii,i}   = Sdivi;
                end
            end
    end

end

end
