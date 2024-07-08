%UI figure for closed loop input (Link input, Connection, Location, Orientation)
%Last modified by Anup Teejo Mathew 02.03.2022
function UIClosedLoop(Tr,iCL)

close all

Tr.plotq0;

fig      = uifigure('Position',[100 100 400 250]);
fig.Name = ['Closed loop joint ',num2str(iCL)];

uilabel(fig,'Position',[25 230 350 25],'Text','Tip of Link A will be connected to tip of Link B');
uilabel(fig,'Position',[25 215 350 25],'Text','(Joint will be defined wrt to the body frame of Link A)');

iList = cell(Tr.N+1,1);
for i=1:Tr.N+1
    iList{i} = num2str(i-1);
end

uilabel(fig,'Position',[25 170 350 25],'Text','Link number corresponding to Link A:');
LinkA = uidropdown(fig,'Position',[230 170 100 22],'Items',iList,'Value','1');
uilabel(fig,'Position',[25 130 350 25],'Text','Link number corresponding to Link B:');
LinkB = uidropdown(fig,'Position',[230 130 100 22],'Items',iList);
uilabel(fig,'Position',[25 80 350 25],'Text','Closed loop joint type:');
CLjAB = uidropdown(fig,'Position',[150 80 100 22],'Items',{'Revolute','Prismatic','Helical','Cylindrical','Planar','Spherical','Fixed'});

% Create push buttons
donebtn = uibutton(fig,'push','Position',[175, 40, 50, 25],'Text','Done',...
                       'ButtonPushedFcn', @(btn,event) DoneButtonPushed(btn,Tr,iCL,LinkA,LinkB,CLjAB));
uiwait(fig)

end

function DoneButtonPushed(btn,Tr,iCL,LinkA,LinkB,CLjAB)
    
    if isequal(LinkA.Value,LinkB.Value)
        uiwait(msgbox('Link A and B are same','Error','error'))
    else
        Tr.nCLj      = Tr.nCLj+1;
        Tr.iACL(iCL) = str2num(LinkA.Value);
        Tr.iCLB(iCL) = str2num(LinkB.Value);
        switch CLjAB.Value
            case 'Revolute'
                L.jointtype='R';
            case 'Prismatic'
                L.jointtype='P';
            case 'Helical'
                L.jointtype='H';
            case 'Cylindrical'
                L.jointtype='C';
            case 'Planar'
                L.jointtype='A';
            case 'Spherical'
                L.jointtype='S';
            case 'Fixed'
                L.jointtype='N';
        end
        
        iA = str2num(LinkA.Value);
        iB = str2num(LinkB.Value);
        g  = Tr.FwdKinematics(zeros(Tr.ndof,1));
        
        i_sigA = 0;
        for i=1:iA
            i_sigA = i_sigA+1;
            if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
                i_sigA = i_sigA+1;
            end
            for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
%                 i_sigA = i_sigA+Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};
                i_sigA = i_sigA+Tr.CVTwists{i}(j+1).nip;
            end
        end
        
        i_sigB = 0;
        for i=1:iB
            i_sigB = i_sigB+1;
            if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
                i_sigB = i_sigB+1;
            end
            for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
                i_sigB = i_sigB++Tr.CVTwists{i}(j+1).nip;
            end
        end
           
        if iA==0
            gBtip = g((i_sigB-1)*4+1:i_sigB*4,:);
            if Tr.VLinks(Tr.LinkIndex(iB)).linktype=='r'
                gBtip = gBtip*Tr.VLinks(Tr.LinkIndex(iB)).gf;
            else
                gBtip = gBtip*Tr.VLinks(Tr.LinkIndex(iB)).gf{end};
            end
            
            gACLj = [eye(3),gBtip(1:3,4);[0 0 0 1]]; %only position
            gBCLj = ginv(gBtip)*gACLj;
        elseif iB==0
            gACLj = eye(4);
            gAtip = g((i_sigA-1)*4+1:i_sigA*4,:);
            if Tr.VLinks(Tr.LinkIndex(iA)).linktype=='r'
                gAtip = gAtip*Tr.VLinks(Tr.LinkIndex(iA)).gf;
            else
                gAtip = gAtip*Tr.VLinks(Tr.LinkIndex(iA)).gf{end};
            end
            gBCLj = gAtip;
        else
            gACLj = eye(4);
            gAtip = g((i_sigA-1)*4+1:i_sigA*4,:);
            if Tr.VLinks(Tr.LinkIndex(iA)).linktype=='r'
                gAtip = gAtip*Tr.VLinks(Tr.LinkIndex(iA)).gf;
            else
                gAtip = gAtip*Tr.VLinks(Tr.LinkIndex(iA)).gf{end};
            end
            gBtip = g((i_sigB-1)*4+1:i_sigB*4,:);
            if Tr.VLinks(Tr.LinkIndex(iB)).linktype=='r'
                gBtip = gBtip*Tr.VLinks(Tr.LinkIndex(iB)).gf;
            else
                gBtip = gBtip*Tr.VLinks(Tr.LinkIndex(iB)).gf{end};
            end
            gBCLj = ginv(gBtip)*gAtip;
        end
        
        Tr.gACLj{iCL} = gACLj;
        Tr.gBCLj{iCL} = gBCLj;
        
        %         Tr.VTwistsCLj(iCL) = jointtwist_pre(L,iCL);
%         Tr.plotq0([Tr.iACL(iCL) Tr.iCLB(iCL)],[],iCL)
%         Tr.VTwistsCLj(iCL) = jointtwist(L,iCL);
%         TwistAB = Tr.VTwistsCLj(iCL);
        TwistAB = jointtwist(L,iCL);
        
        save('Temp_LinkageClosedJoint','iA','iB','TwistAB','gACLj','gBCLj')
        close all
        closereq(); 
    end
end
