%UI figure for link input (Link input, Connection, Location, Orientation)
%Last modified by Anup Teejo Mathew 02.03.2022
function UILinkInput(Tr,LinkNames,i_lnk)

close all

Tr.plotq0;

fig      = uifigure('Position',[100 100 400 325]);
fig.Name = sprintf('Link %d input',i_lnk);

uilabel(fig,'Position',[25 275 350 25],'Text','Link name:');
LinkName = uidropdown(fig,'Position',[150 275 100 22],'Items',LinkNames);

bg  = uibuttongroup(fig,'Position',[25 175 200 30]);
rb1 = uiradiobutton(bg,'Position',[5 5 95 25]);
rb2 = uiradiobutton(bg,'Position',[105 5 95 25]);

rb1.Text = 'Location';
rb2.Text = 'Orientation';

text1 = sprintf('Link %d location (in m) (orientation (cartesian plane angles in rads))',i_lnk);
uilabel(fig,'Position',[25 130 350 50],'Text',text1);
l1 = uilabel(fig,'Position',[25 115 350 50],'Text','(With respect to the body frame at the tip of previous Link)');

uilabel(fig,'Position',[30 100 35 25],'Text','x (qx):');
x = uieditfield(fig,'Position',[75 100 50 25],'Value','0');

uilabel(fig,'Position',[130 100 35 25],'Text','y (qy):');
y = uieditfield(fig,'Position',[175 100 50 25],'Value','0');

uilabel(fig,'Position',[230 100 35 25],'Text','z (qz):');
z = uieditfield(fig,'Position',[275 100 50 25],'Value','0');

if Tr.N>0
    
    uilabel(fig,'Position',[25 225 350 25],'Text','Link number to which this link is connected:');
    uilabel(fig,'Position',[25 210 350 25],'Text','(Previous Link):');
    iList = cell(Tr.N,1);
    for i=1:Tr.N+1
        iList{i} = num2str(i-1);
    end
    PreviousLink = uidropdown(fig,'Position',[275 225 100 22],'Items',iList,'Value',num2str(Tr.N),...
                                  'ValueChangedFcn',@(PreviousLink,event) numberChanged(PreviousLink,l1));
else
    PreviousLink.Value='0';
end

global g_ini_i
g_ini_i=eye(4);

% Create push buttons
applybtn = uibutton(fig,'push','Position',[75, 50, 50, 25],'Text','Apply',...
                        'ButtonPushedFcn', @(btn,event) ApplyButtonPushed(btn,Tr,LinkNames,LinkName,PreviousLink,rb1,x,y,z,l1));
donebtn = uibutton(fig,'push','Position',[175, 50, 50, 25],'Text','Done',...
                       'ButtonPushedFcn', @(btn,event) DoneButtonPushed(btn,LinkNames,LinkName,PreviousLink,rb1,x,y,z));
resetbtn = uibutton(fig,'push','Position',[275, 50, 50, 25],'Text','Reset',...
                        'ButtonPushedFcn', @(btn,event) ResetButtonPushed(btn,Tr,LinkNames,LinkName,PreviousLink,x,y,z));
uiwait(fig)

end

% Create the function for the ButtonPushedFcn callback
function ApplyButtonPushed(btn,Tr,LinkNames,LinkName,PreviousLink,rb1,x,y,z,l1)

global g_ini_i
if rb1.Value
    xi      = [0;0;0;str2num(x.Value);str2num(y.Value);str2num(z.Value)];
    g_ini_i = g_ini_i*variable_expmap_g(xi);
    x.Value = '0'; y.Value = '0'; z.Value = '0';
else
    xi      = [str2num(x.Value);str2num(y.Value);str2num(z.Value);0;0;0];
    g_ini_i = g_ini_i*variable_expmap_g(xi);
    x.Value = '0'; y.Value = '0'; z.Value = '0';
end

Tr.N        = Tr.N+1;
LinkIndex_i = find(strcmp(LinkNames,LinkName.Value));
iLpre_i     = str2num(PreviousLink.Value);


Tr.LinkIndex = [Tr.LinkIndex;LinkIndex_i];
Tr.iLpre     = [Tr.iLpre;iLpre_i];
Tr.g_ini     = [Tr.g_ini;g_ini_i];

Lscale = 0;
for i=1:Tr.N
    Lscale = Lscale+Tr.VLinks(Tr.LinkIndex(i)).L;
end
Tr.PlotParameters.Lscale = Lscale;

ndof    = Tr.ndof;
VTwists = SorosimTwist.empty(Tr.VLinks(LinkIndex_i).npie,0);

T          = jointtwist_pre(Tr.VLinks(LinkIndex_i),[]); %for each joint and rigid link
VTwists(1) = T;
ndof       = ndof+T.dof;

for j=1:Tr.VLinks(LinkIndex_i).npie-1 %for each of the soft link divisions

    T            = SorosimTwist;
    VTwists(j+1) = T;

end

Tr.CVTwists = [Tr.CVTwists;{VTwists}];
Tr.ndof     = ndof;

close(1)
Tr.plotq0;
l1.Text = '(With respect to the body frame at the base of the Link)';

end

function DoneButtonPushed(btn,LinkNames,LinkName,PreviousLink,rb1,x,y,z)

global g_ini_i
if rb1.Value
    xi      = [0;0;0;str2num(x.Value);str2num(y.Value);str2num(z.Value)];
    g_ini_i = g_ini_i*variable_expmap_g(xi);
else
    xi      = [str2num(x.Value);str2num(y.Value);str2num(z.Value);0;0;0];
    g_ini_i = g_ini_i*variable_expmap_g(xi);
end

LinkIndex_i = find(strcmp(LinkNames,LinkName.Value));
iLpre_i     = str2num(PreviousLink.Value);

save('Temp_LinkageAssembly','g_ini_i','LinkIndex_i','iLpre_i')
close all
closereq();
end

function ResetButtonPushed(btn,Tr,LinkNames,LinkName,PreviousLink,x,y,z)

global g_ini_i
g_ini_i = eye(4);
x.Value = '0'; y.Value = '0'; z.Value = '0';

Tr.N        = Tr.N+1;
LinkIndex_i = find(strcmp(LinkNames,LinkName.Value));
iLpre_i     = str2num(PreviousLink.Value);

Tr.LinkIndex = [Tr.LinkIndex;LinkIndex_i];
Tr.iLpre     = [Tr.iLpre;iLpre_i];
Tr.g_ini     = [Tr.g_ini;g_ini_i];

Lscale = 0;
for i=1:Tr.N
    Lscale = Lscale+Tr.VLinks(Tr.LinkIndex(i)).L;
end
Tr.PlotParameters.Lscale = Lscale;

ndof    = Tr.ndof;
VTwists = SorosimTwist.empty(Tr.VLinks(LinkIndex_i).npie,0);

T          = jointtwist_pre(Tr.VLinks(LinkIndex_i),[]); %for each joint and rigid link
VTwists(1) = T;
ndof       = ndof+T.dof;

for j=1:Tr.VLinks(LinkIndex_i).npie-1 %for each of the soft link divisions

    T            = Twist;
    VTwists(j+1) = T;

end

Tr.CVTwists = [Tr.CVTwists;{VTwists}];
Tr.ndof     = ndof;
Tr.plotq0;
end

function numberChanged(PreviousLink,l1)
l1.Text = '(With respect to the body frame at the tip of previous Link)';
global g_ini_i
g_ini_i = eye(4);
end
