%Function to input the value of joint wrenches/coordinates
%Last modified by Anup Teejo Mathew 02.03.2022
function uq = InputJointUQ0(Tr)

n_jact             = Tr.n_jact;
i_jact             = Tr.i_jact;
uq                 = zeros(n_jact,1);
WrenchControlled   = Tr.WrenchControlled;

i_qj      = 1;
VLinks    = Tr.VLinks;
LinkIndex = Tr.LinkIndex;

for i=i_jact
    if VLinks(LinkIndex(i)).jointtype=='R'

        if WrenchControlled(i_qj)
            prompt   = ['Enter the value of torque at the revolute joint of link ',num2str(i),' (Nm):'];
            answer   = input(prompt, 's');
            uq(i_qj) = str2num(answer);
        else
            prompt   = ['Enter the value of angle at the revolute joint of link ',num2str(i),' (rad):'];
            answer   = input(prompt, 's');
            uq(i_qj) = str2num(answer);
        end
        i_qj=i_qj+1;

    elseif VLinks(LinkIndex(i)).jointtype=='P'

        if WrenchControlled(i_qj)
            prompt   = ['Enter the value of force at the prismatic joint of link ',num2str(i),' (N):'];
            answer   = input(prompt, 's');
            uq(i_qj) = str2num(answer);
        else
            prompt   = ['Enter the value of displacement at the prismatic joint of link ',num2str(i),' (m):'];
            answer   = input(prompt, 's');
            uq(i_qj) = str2num(answer);
        end
        i_qj=i_qj+1;

    elseif VLinks(LinkIndex(i)).jointtype=='H'

        if WrenchControlled(i_qj)
            prompt={['Enter the value of torque at the helical joint of link ',num2str(i),' (Nm):']...
                ,['Enter the value of force at the helical joint of link ',num2str(i),' (N):']};
            dlgtitle         = 'Wrench controlled helical joint';
            definput         = {'0.1','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            f=1;
            for ii=1:i-1
                f = f+Tr.VLinks(ii).npie;
            end
            B        = Tr.Vtwists(f).B;
            pitch    = sum(B(4:6));
            uq(i_qj) = str2num(answer{1})+pitch*str2num(answer{2});
        else
            prompt   =['Enter the value of angle at the helical joint of link ',num2str(i),' as a function of t (m):'];
            answer   = input(prompt, 's');
            uq(i_qj) = str2num(answer);
        end
        i_qj=i_qj+1;
    elseif VLinks(LinkIndex(i)).jointtype=='U'

        if WrenchControlled(i_qj)
            prompt={['Enter the value of torque 1 at the universal joint of link ',num2str(i),' (Nm):']...
                ,['Enter the value of torque 2 at the universal joint of link ',num2str(i),' (Nm):']};
            dlgtitle         = 'Wrench controlled universal joint';
            definput         = {'0.1','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
        else
            prompt={['Enter the value of angle 1 at the universal joint of link ',num2str(i),' (rad):']...
                ,['Enter the value of angle 2 at the universal joint of link ',num2str(i),' (rad):']};
            dlgtitle         = 'Joint coordinate controlled universal joint';
            definput         = {'pi/2','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
        end
        i_qj=i_qj+2;

    elseif VLinks(LinkIndex(i)).jointtype=='C'

        if WrenchControlled(i_qj)
            prompt={['Enter the value of torque at the cylindrical joint of link ',num2str(i),' (Nm):']...
                ,['Enter the value of force at the cylindrical joint of link ',num2str(i),' (N):']};
            dlgtitle         = 'Wrench controlled cylindrical joint';
            definput         = {'0.1','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
        else
            prompt={['Enter the value of angle at the cylindrical joint of link ',num2str(i),' (rad):']...
                ,['Enter the value of displacement at the cylindrical joint of link ',num2str(i),' (m):']};
            dlgtitle         = 'Joint coordinate controlled cylindrical joint';
            definput         = {'pi/2','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
        end
        i_qj=i_qj+2;
    elseif VLinks(LinkIndex(i)).jointtype=='A'

        if WrenchControlled(i_qj)
            prompt={['Enter the value of torque at the planar joint of link ',num2str(i),' (Nm):']...
                ,['Enter the value of force 1 at the planar joint of link ',num2str(i),' (N):']...
                ,['Enter the value of force 2 at the planar joint of link ',num2str(i),' (N):']};
            dlgtitle         = 'Wrench controlled planar joint';
            definput         = {'0.1','0','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
            uq(i_qj+2)  = str2num(answer{3});
        else
            prompt={['Enter the value of angle at the planar joint of link ',num2str(i),' (rad):']...
                ,['Enter the value of displacement 1 at the planar joint of link ',num2str(i),' (m):']...
                ,['Enter the value of displacement 2 at the planar joint of link ',num2str(i),' (m):']};
            dlgtitle = 'Joint coordinate controlled planar joint';

            definput = {'pi/2','0','0'};
            opts.Interpreter = 'tex';
            answer = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
            uq(i_qj+2)  = str2num(answer{3});
        end
        i_qj=i_qj+3;

    elseif VLinks(LinkIndex(i)).jointtype=='S'

        if WrenchControlled(i_qj)
            prompt={['Enter the value of torque (Mx) at the spherical joint of link ',num2str(i),' (Nm):']...
                ,['Enter the value of torque (My) at the spherical joint of link ',num2str(i),' (Nm):']...
                ,['Enter the value of torque (Mz) at the spherical joint of link ',num2str(i),' (Nm):']};
            dlgtitle         = 'Wrench controlled spherical joint';

            definput         = {'0.1','0','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
            uq(i_qj+2)  = str2num(answer{3});
        else
            prompt={['Enter the value of q1 at the spherical joint of link ',num2str(i),' (rad):']...
                ,['Enter the value of q2 at the spherical joint of link ',num2str(i),' (rad):']...
                ,['Enter the value of q3 at the spherical joint of link ',num2str(i),' (rad):']};
            dlgtitle = 'Joint coordinate controlled spherical joint';

            definput = {'pi/2','0','0'};
            opts.Interpreter = 'tex';
            answer = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
            uq(i_qj+2)  = str2num(answer{3});
        end
        i_qj=i_qj+3;

    elseif VLinks(LinkIndex(i)).jointtype=='F'

        if WrenchControlled(i_qj)
            prompt={['Enter the value of torque (Mx) at the free joint of link ',num2str(i),' (Nm):']...
                ,['Enter the value of torque (My) at the free joint of link ',num2str(i),' (Nm):']...
                ,['Enter the value of torque (Mz) at the free joint of link ',num2str(i),' (Nm):']...
                ,['Enter the value of force (Fx) at the free joint of link ',num2str(i),' (N):']...
                ,['Enter the value of force (Fy) at the free joint of link ',num2str(i),' (N):']...
                ,['Enter the value of force (Fz) at the free joint of link ',num2str(i),' (N):']};
            dlgtitle = 'Wrench controlled free joint';

            definput = {'0.1','0','0','0','0','0'};
            opts.Interpreter = 'tex';
            answer = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
            uq(i_qj+2)  = str2num(answer{3});
            uq(i_qj+3)  = str2num(answer{4});
            uq(i_qj+4)  = str2num(answer{5});
            uq(i_qj+5)  = str2num(answer{6});
        else
            prompt={['Enter the value of q1 at the free joint of link ',num2str(i),' (rad):']...
                ,['Enter the value of q2 at the free joint of link ',num2str(i),' (rad):']...
                ,['Enter the value of q3 at the free joint of link ',num2str(i),' (rad):']...
                ,['Enter the value of q4 at the free joint of link ',num2str(i),' (m):']...
                ,['Enter the value of q5 at the free joint of link ',num2str(i),' (m):']...
                ,['Enter the value of q6 at the free joint of link ',num2str(i),' (m):']};
            dlgtitle = 'Joint coordinate controlled free joint';

            definput = {'pi/2','0','0','0','0','0'};
            opts.Interpreter = 'tex';
            answer = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uq(i_qj)    = str2num(answer{1});
            uq(i_qj+1)  = str2num(answer{2});
            uq(i_qj+2)  = str2num(answer{3});
            uq(i_qj+3)  = str2num(answer{4});
            uq(i_qj+4)  = str2num(answer{5});
            uq(i_qj+5)  = str2num(answer{6});
        end
        i_qj=i_qj+6;

    else
        error('Check joint actuation');
    end
end
end

