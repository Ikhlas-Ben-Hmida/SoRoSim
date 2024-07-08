%Function to Assign Twist class (basis, dof, etc.) to joints
%Last modified by Anup Teejo Mathew 02.03.2022

function T = jointtwist(L,i) 

    switch L.jointtype
        case 'R' %if revolute
            
            list  = {'x', 'y', 'z'};
            J_axis = listdlg('PromptString',{['Revolute joint of link ',num2str(i)],'Axis of rotation (local):'},...
                             'SelectionMode','single','ListSize',[160 160],'ListString',list,'InitialValue',3);
            if J_axis==1
                B = [1 0 0 0 0 0]';
            elseif J_axis==2
                B = [0 1 0 0 0 0]';
            else
                B = [0 0 1 0 0 0]';
            end

            T        = SorosimTwist(B);
            T.dof    = 1;

        case 'P' %if prismatic
            
            list  = {'x', 'y', 'z'};
            J_axis = listdlg('PromptString',{['Prismatic joint of link ',num2str(i)],'Axis of translation (local):'},...
                             'SelectionMode','single','ListSize',[160 160],'ListString',list,'InitialValue',3);
            if J_axis==1
                B = [0 0 0 1 0 0]';
            elseif J_axis==2
                B = [0 0 0 0 1 0]';
            else
                B = [0 0 0 0 0 1]';
            end

            T        = SorosimTwist(B);
            T.dof    = 1;

        case 'H' %if helical
            
            badans = 1;
            
            while badans
                badans=0;
            
                prompt           = {'Joint pitch (m/rad):','Axis of motion (x, y or z) (local):'};
                dlgtitle         = ['Helical joint of link ',num2str(i)];
                definput         = {'0.1','x'};
                opts.Interpreter = 'tex';
                opts.WindowStyle = 'Normal';
                H_ans            = inputdlg(prompt,dlgtitle,[1 75],definput,opts);

                p      = str2num(H_ans{1});
                J_axis = H_ans{2};
                if J_axis=='x'
                    B = [1 0 0 p 0 0]';
                elseif J_axis== 'y'
                    B = [0 1 0 0 p 0]';
                elseif J_axis== 'z'
                    B = [0 0 1 0 0 p]';
                else
                    uiwait(msgbox('Incorrect entry of axis','Error','error'));
                    badans=1;
                end
            end

            T        = SorosimTwist(B);
            T.dof    = 1;
            

        case 'C' %if cylinderical
            
            list  = {'x', 'y', 'z'};
            J_axis = listdlg('PromptString',{['Cylinderical joint of link ',num2str(i)],'Axis of motion (local):'},...
                             'SelectionMode','single','ListSize',[160 160],'ListString',list,'InitialValue',2);
                
            if J_axis==1
                B = [1 0 0 0 0 0;0 0 0 1 0 0]'; %correct this later
            elseif J_axis== 2
                B = [0 1 0 0 0 0;0 0 0 0 1 0]';
            else
                B = [0 0 1 0 0 0;0 0 0 0 0 1]';
            end

            T        = SorosimTwist(B);
            T.dof    = 2;

        case 'A' %if planar

            
            list  = {'xy', 'yz', 'xz'};
            J_axis = listdlg('PromptString',{['Planar joint of link ',num2str(i)],'Plane of motion (local):'},...
                             'SelectionMode','single','ListSize',[160 160],'ListString',list,'InitialValue',1);

            if J_axis==1
                B = [0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0]';
            elseif J_axis==2
                B = [1 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 1]';
            else
                B = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]';
            end

            T     = SorosimTwist(B);
            T.dof = 3;

        case 'S' %if spherical
            
            B = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0]';

            T        = SorosimTwist(B);
            T.dof    = 3;

        case 'F' %if free motion

            B=eye(6);

            T        = SorosimTwist(B);
            T.dof    = 6;
            
        case 'N'
            
            B     = [];
            T     = SorosimTwist(B);
            T.dof = 0;
            
        otherwise
            error(['Incorrect joint type. Check the joint type of the link ',num2str(i)])
    end
    T.xi_star=zeros(6,1);
    
end

