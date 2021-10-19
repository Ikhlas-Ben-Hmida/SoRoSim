%Function to Assign Twist class (basis, dof, etc.) to joints
%Last modified by Anup Teejo Mathew 23/05/2021

function T = jointtwist(L,i) 

    switch L.jointtype
        case 'R' %if revolute
            
            prompt           = {'Axis of rotation (x,y or z):'};
            dlgtitle         = ['Revolute joint of link ',num2str(i)];
            definput         = {'z'};
            opts.Interpreter = 'tex';
            opts.WindowStyle = 'Normal';
            R_ans            = inputdlg(prompt,dlgtitle,[1 75],definput,opts);
            
            J_axis = R_ans{1};
            if J_axis=='x'
                B = [1 0 0 0 0 0]';
            elseif J_axis=='y'
                B = [0 1 0 0 0 0]';
            else
                B = [0 0 1 0 0 0]';
            end

            T        = Twist(i,1,[],0,B);
            T.dof    = 1;

        case 'P' %if prismatic
            
            prompt           = {'Axis of translation (x,y or z):'};
            dlgtitle         = ['Prismatic joint of link ',num2str(i)];
            definput         = {'x'};
            opts.Interpreter = 'tex';
            opts.WindowStyle = 'Normal';
            P_ans            = inputdlg(prompt,dlgtitle,[1 75],definput,opts);
            
            J_axis = P_ans{1};
            if J_axis=='x'
                B = [0 0 0 1 0 0]';
            elseif J_axis== 'y'
                B = [0 0 0 0 1 0]';
            else
                B = [0 0 0 0 0 1]';
            end

            T        = Twist(i,1,[],0,B);
            T.dof    = 1;

        case 'H' %if helical
            
            prompt           = {'Joint pitch (m/rad):','Axis of motion:'};
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
            else
                B = [0 0 1 0 0 p]';
            end

            T        = Twist(i,1,[],0,B);
            T.dof    = 1;
            
        case 'U' %if universal
            
            prompt           = {'Constrained axis of rotation:'};
            dlgtitle         = ['Universal joint of link ',num2str(i)];
            definput         = {'x'};
            opts.Interpreter = 'tex';
            opts.WindowStyle = 'Normal';
            C_ans            = inputdlg(prompt,dlgtitle,[1 75],definput,opts);

            J_axis1 = C_ans{1};
            if J_axis1=='x'
                B = [0 1 0 0 0 0;0 0 1 0 0 0]'; 
            elseif J_axis1== 'y'
                B = [1 0 0 0 0 0;0 0 1 0 0 0]'; 
            else
                B = [1 0 0 0 0 0;0 1 0 0 0 0]'; 
            end

            T        = Twist(i,1,[],0,B);
            T.dof    = 2;

        case 'C' %if cylinderical
            
            prompt           = {'Axis of motion:'};
            dlgtitle         = ['Cylinderical joint of link ',num2str(i)];
            definput         = {'y'};
            opts.Interpreter = 'tex';
            opts.WindowStyle = 'Normal';
            C_ans            = inputdlg(prompt,dlgtitle,[1 75],definput,opts);

            J_axis1 = C_ans{1};
            if J_axis1=='x'
                B = [1 0 0 0 0 0;0 0 0 1 0 0]'; %correct this later
            elseif J_axis1== 'y'
                B = [0 1 0 0 0 0;0 0 0 0 1 0]';
            else
                B = [0 0 1 0 0 0;0 0 0 0 0 1]';
            end

            T        = Twist(i,1,[],0,B);
            T.dof    = 2;

        case 'A' %if planar
            
            prompt           = {'Plane of motion: (xy, yz or xz)'};
            dlgtitle         = ['Planar joint of link ',num2str(i)];
            definput         = {'xy'};
            opts.Interpreter = 'tex';
            opts.WindowStyle = 'Normal';
            A_ans            = inputdlg(prompt,dlgtitle,[1 75],definput,opts);

            J_axis = A_ans{1};
            if all(J_axis=='xy')
                B = [0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0]';
            elseif all(J_axis== 'yz')
                B = [1 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 1]';
            else
                B = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]';
            end

            T=Twist(i,1,[],0,B);
            T.dof=3;

        case 'S' %if spherical
            
            B = [1 0 0 0 0 0;0 1 0 0 0 0;0 0 1 0 0 0]';

            T        = Twist(i,1,[],0,B);
            T.dof    = 3;

        case 'F' %if free motion

            B=eye(6);

            T        = Twist(i,1,[],0,B);
            T.dof    = 6;
            
        case 'N'
            
            B     = [];
            T     = Twist(i,1,[],0,B);
            T.dof = 0;
            
        otherwise
            error(['Incorrect joint type. Check the joint type of the link ',num2str(i)])
    end
    T.xi_star=zeros(6,1);
    
end

