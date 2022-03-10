%Function that allows the user to set actuation cable parameters
%(24.05.2021)
function [n_sact,dc,dcp,Sdiv,Ediv,Inside] = CableActuation(Tr)

prompt           = {'Number of actuators:'};
dlgtitle         = 'Actuation Parameters';
definput         = {'1'};
opts.Interpreter = 'tex';
ans_act          = inputdlg(prompt,dlgtitle,[1 50],definput,opts);
n_sact           = str2num(ans_act{1});
N                = Tr.N;

dc               = cell(n_sact,N);
dcp              = cell(n_sact,N);
Sdiv             = cell(n_sact,N);
Ediv             = cell(n_sact,N);
Inside           = cell(1,n_sact);

Tr.n_sact        = n_sact;

%Set cable parameters
for ii = 1:n_sact
    
    if any(Tr.CS=='R')
        quest  = ['Is the actutator (cable) ',num2str(ii),' fully inside the linkage?'];
        Answer = questdlg(quest,'Inside?','Yes','No','Yes');
    else
        Answer = 'Yes';
    end

    switch Answer
        case 'Yes'
            Inside{ii} = true;
            for i=1:N
                if Tr.VLinks(Tr.LinkIndex(i)).linktype=='s'
                    cableactuationGUI(Tr,ii,i)
                    load('cableactuation.mat','dci','dcpi','div_starti','div_endi')
                    dc{ii,i}     = dci;
                    dcp{ii,i}    = dcpi;
                    Ediv{ii,i}   = div_endi;
                    Sdiv{ii,i}   = div_starti;
                end
            end
        case 'No'
            Inside{ii} = false;
            for i=1:N
                if Tr.VLinks(Tr.LinkIndex(i)).linktype=='s'
                    cableactuation2GUI(Tr,ii,i)
                    load('cableactuation.mat','dci','div_starti','div_endi')
                    dc{ii,i}     = dci;
                    Ediv{ii,i}   = div_endi;
                    Sdiv{ii,i}   = div_starti;
                end
            end
    end

end

end
