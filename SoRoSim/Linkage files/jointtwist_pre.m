%Function to Pre assign Twist class (basis, dof, etc.) to joints to be
%changed during Linkage creation
%Last modified by Anup Teejo Mathew 30/05/2021

function T = jointtwist_pre(L,i) 

    switch L.jointtype
        case 'R' %if revolute

            J_axis = 'z';
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
            
            J_axis = 'x';
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
            
            p      = 0.1;
            J_axis = 'x';
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

            J_axis1 = 'x';
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
            
            J_axis1 = 'y';
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

            J_axis = 'xy';
            if all(J_axis=='xy')
                B = [0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 0 1 0]';
            elseif all(J_axis== 'yz')
                B = [1 0 0 0 0 0;0 0 0 0 1 0;0 0 0 0 0 1]';
            else
                B = [0 1 0 0 0 0;0 0 0 1 0 0;0 0 0 0 0 1]';
            end

            T        = Twist(i,1,[],0,B);
            T.dof    = 3;

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

