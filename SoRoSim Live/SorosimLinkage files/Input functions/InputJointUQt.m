%Function to input the value of joint wrenches/coordinates as function of time
%Last modified by Anup Teejo Mathew 02.03.2022
function uqt = InputJointUQt(Tr)

n_jact             = Tr.n_jact;
i_jact             = Tr.i_jact;
uqt                = cell(n_jact,1);
WrenchControlled   = Tr.WrenchControlled;

i_qj      = 1;
VLinks    = Tr.VLinks;
LinkIndex = Tr.LinkIndex;

for i=i_jact
    
    if VLinks(LinkIndex(i)).jointtype=='R'
        
        if WrenchControlled(i_qj)
            prompt    = ['Enter the value of torque at the revolute joint of link ',num2str(i),' as a function of t (Nm):'...
                         '\n[Examples: 0.1*t, 0.1*sin(2*pi*t)]: '];
            funstr    = input(prompt, 's');
            uqt{i_qj} = str2func( ['@(t) ' funstr ] );
        else
            prompt = ['Enter the value of angle at the revolute joint of link ',num2str(i),' as a function of t (rad):'...
                      '\n[Examples: 2*pi*t, pi/2*sin(2*pi*t)]: '];
            funstr = input(prompt, 's');
            q      = str2sym(funstr);
            qd     = diff(q);
            qdd    = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))] );
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj}={q;qd;qdd};
        end
        i_qj=i_qj+1;
        
    elseif VLinks(LinkIndex(i)).jointtype=='P'
        
        if WrenchControlled(i_qj)
            prompt    = ['Enter the value of force at the prismatic joint of link ',num2str(i),' as a function of t (N):'...
                         '\n[Examples: 2*t, 2*sin(2*pi*t)]: '];
            funstr    = input(prompt, 's');
            uqt{i_qj} = str2func( ['@(t) ' funstr ] );
        else
            prompt = ['Enter the value of displacement at the prismatic joint of link ',num2str(i),' as a function of t (m):'...
                      '\n[Examples: 0.1*t, 0.2*sin(2*pi*t)]: '];
            funstr = input(prompt, 's');
            
            q    = str2sym(funstr);
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj}={q;qd;qdd};
        end
        i_qj=i_qj+1;
        
    elseif VLinks(LinkIndex(i)).jointtype=='H'
        
        if WrenchControlled(i_qj)
            prompt = {['Enter the value of torque at the helical joint of link ',num2str(i),' as a function of t (Nm):']...
                   ,['Enter the value of force at the helical joint of link ',num2str(i),' as a function of t (N):']};
            dlgtitle         = 'Wrench controlled helical joint';
            definput         = {'0.1*t','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            f=1;
            for ii=1:i-1
                f=f+Tr.VLinks(ii).npie;
            end
            B         = Tr.Vtwists(f).B;
            pitch     = sum(B(4:6));
            uH        = [answer{1},'+',num2str(pitch),'*(',answer{2},')'];
            uqt{i_qj} = str2func( ['@(t) ' uH ] );
        else
            prompt = ['Enter the value of angle at the helical joint of link ',num2str(i),' as a function of t (m):'];
            funstr = input(prompt, 's');
            
            q    = str2sym(funstr);
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj}={q;qd;qdd};
        end
        i_qj=i_qj+1;
        
    elseif VLinks(LinkIndex(i)).jointtype=='U'
        
         if WrenchControlled(i_qj)
            prompt = {['Enter the value of torque 1 at the universal joint of link ',num2str(i),' as a function of t (Nm):']...
                   ,['Enter the value of torque 2 at the universal joint of link ',num2str(i),' as a function of t (Nm):']};
            dlgtitle         = 'Wrench controlled universal joint';
            definput         = {'0.1*t','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uqt{i_qj}    = str2func( ['@(t) ' answer{1} ] );
            uqt{i_qj+1}  = str2func( ['@(t) ' answer{2} ] );
        else
            prompt={['Enter the value of angle 1 at the universal joint of link ',num2str(i),' as a function of t (rad):']...
                   ,['Enter the value of angle 2 at the universal joint of link ',num2str(i),' as a function of t (rad):']};
            dlgtitle         = 'Joint coordinate controlled universal joint';
            definput         = {'2*pi*t','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            
            %1
            
            q    = str2sym(answer{1});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj}={q;qd;qdd};
            
            %2
            
            q    = str2sym(answer{2});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+1}={q;qd;qdd};
        end
        i_qj=i_qj+2;
        
    elseif VLinks(LinkIndex(i)).jointtype=='C'
        
        if WrenchControlled(i_qj)
            prompt = {['Enter the value of torque at the cylindrical joint of link ',num2str(i),' as a function of t (Nm):']...
                   ,['Enter the value of force at the cylindrical joint of link ',num2str(i),' as a function of t (N):']};
            dlgtitle         = 'Wrench controlled cylindrical joint';
            definput         = {'0.1*t','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            uqt{i_qj}    = str2func( ['@(t) ' answer{1} ] );
            uqt{i_qj+1}  = str2func( ['@(t) ' answer{2} ] );
        else
            prompt={['Enter the value of angle at the cylindrical joint of link ',num2str(i),' as a function of t (rad):']...
                   ,['Enter the value of displacement at the cylindrical joint of link ',num2str(i),' as a function of t (m):']};
            dlgtitle         = 'Joint coordinate controlled cylindrical joint';
            definput         = {'2*pi*t','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            
            %1
            
            q    = str2sym(answer{1});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj}={q;qd;qdd};
            
            %2
            
            q    = str2sym(answer{2});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+1}={q;qd;qdd};
        end
        i_qj=i_qj+2;
    elseif VLinks(LinkIndex(i)).jointtype=='A'
        
        if WrenchControlled(i_qj)
            prompt = {['Enter the value of torque at the planar joint of link ',num2str(i),' as a function of t (Nm):']...
                   ,['Enter the value of force 1 at the planar joint of link ',num2str(i),' as a function of t (N):']...
                   ,['Enter the value of force 2 at the planar joint of link ',num2str(i),' as a function of t (N):']};
            dlgtitle         = 'Wrench controlled planar joint';
            definput         = {'0.1*t','0','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            
            uqt{i_qj}    = str2func( ['@(t) ' answer{1} ] );
            uqt{i_qj+1}  = str2func( ['@(t) ' answer{2} ] );
            uqt{i_qj+2}  = str2func( ['@(t) ' answer{3} ] );
            
        else
            prompt = {['Enter the value of angle at the planar joint of link ',num2str(i),' as a function of t (rad):']...
                   ,['Enter the value of displacement 1 at the planar joint of link ',num2str(i),' as a function of t (m):']...
                   ,['Enter the value of displacement 2 at the planar joint of link ',num2str(i),' as a function of t (m):']};
            dlgtitle         = 'Joint coordinate controlled planar joint';
            definput         = {'2*pi*t','0','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            
            %1
            
            q    = str2sym(answer{1});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj}={q;qd;qdd};
            
            %2
            
            q    = str2sym(answer{2});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+1}={q;qd;qdd};
            
            %3
            
            q    = str2sym(answer{3});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+2}={q;qd;qdd};
        end
        i_qj=i_qj+3;
              
    elseif VLinks(LinkIndex(i)).jointtype=='S'
        
        if WrenchControlled(i_qj)
            prompt = {['Enter the value of torque (Mx) at the spherical joint of link ',num2str(i),' as a function of t (Nm):']...
                   ,['Enter the value of torque (My) at the spherical joint of link ',num2str(i),' as a function of t (Nm):']...
                   ,['Enter the value of torque (Mz) at the spherical joint of link ',num2str(i),' as a function of t (Nm):']};
            dlgtitle         = 'Wrench controlled spherical joint';
            definput         = {'0.1*t','0','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            
            uqt{i_qj}    = str2func( ['@(t) ' answer{1} ] );
            uqt{i_qj+1}  = str2func( ['@(t) ' answer{2} ] );
            uqt{i_qj+2}  = str2func( ['@(t) ' answer{3} ] );
        else
            prompt = {['Enter the value of q1 at the spherical joint of link ',num2str(i),' as a function of t (rad):']...
                   ,['Enter the value of q2 at the spherical joint of link ',num2str(i),' as a function of t (rad):']...
                   ,['Enter the value of q3 at the spherical joint of link ',num2str(i),' as a function of t (rad):']};
            dlgtitle         = 'Joint coordinate controlled spherical joint';
            definput         = {'2*pi*t','0','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            
            %1
            
            q    = str2sym(answer{1});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj}={q;qd;qdd};
            
            %2
            
            q    = str2sym(answer{2});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+1}={q;qd;qdd};
            
            %3
            
            q    = str2sym(answer{3});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+2}={q;qd;qdd};
        end
        i_qj=i_qj+3;
        
    elseif VLinks(LinkIndex(i)).jointtype=='F'
        
        if WrenchControlled(i_qj)
            prompt = {['Enter the value of torque (Mx) at the free joint of link ',num2str(i),' as a function of t (Nm):']...
                   ,['Enter the value of torque (My) at the free joint of link ',num2str(i),' as a function of t (Nm):']...
                   ,['Enter the value of torque (Mz) at the free joint of link ',num2str(i),' as a function of t (Nm):']...
                   ,['Enter the value of force (Fx) at the free joint of link ',num2str(i),' as a function of t (N):']...
                   ,['Enter the value of force (Fy) at the free joint of link ',num2str(i),' as a function of t (N):']...
                   ,['Enter the value of force (Fz) at the free joint of link ',num2str(i),' as a function of t (N):']};
            dlgtitle         = 'Wrench controlled free joint';
            definput         = {'0.1*t','0','0','0','0','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            
            uqt{i_qj}    = str2func( ['@(t) ' answer{1} ] );
            uqt{i_qj+1}  = str2func( ['@(t) ' answer{2} ] );
            uqt{i_qj+2}  = str2func( ['@(t) ' answer{3} ] );
            uqt{i_qj+3}  = str2func( ['@(t) ' answer{4} ] );
            uqt{i_qj+4}  = str2func( ['@(t) ' answer{5} ] );
            uqt{i_qj+5}  = str2func( ['@(t) ' answer{6} ] );
        else
            prompt = {['Enter the value of q1 at the free joint of link ',num2str(i),' as a function of t (rad):']...
                   ,['Enter the value of q2 at the free joint of link ',num2str(i),' as a function of t (rad):']...
                   ,['Enter the value of q3 at the free joint of link ',num2str(i),' as a function of t (rad):']...
                   ,['Enter the value of q4 at the free joint of link ',num2str(i),' as a function of t (m):']...
                   ,['Enter the value of q5 at the free joint of link ',num2str(i),' as a function of t (m):']...
                   ,['Enter the value of q6 at the free joint of link ',num2str(i),' as a function of t (m):']};
            dlgtitle         = 'Joint coordinate controlled free joint';
            definput         = {'2*pi*t','0','0','0','0','0'};
            opts.Interpreter = 'tex';
            answer           = inputdlg(prompt,dlgtitle,[1 100],definput,opts);
            
            %1
            
            q    = str2sym(answer{1});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj}={q;qd;qdd};
            
            %2
            
            q    = str2sym(answer{2});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+1}={q;qd;qdd};
            
            %3
            
            q    = str2sym(answer{3});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+2}={q;qd;qdd};
            
            %4
            
            q    = str2sym(answer{4});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+3}={q;qd;qdd};
            
            %5
            
            q    = str2sym(answer{5});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+4}={q;qd;qdd};
            
            %6
            
            q    = str2sym(answer{6});
            qd   = diff(q);
            qdd  = diff(qd);
            
            syms t
            if has(q,t)
                q=matlabFunction(q);
            else
                q=str2func(['@(t) ' num2str(double(q))]);
            end

            if has(qd,t)
                qd=matlabFunction(qd);
            else
                qd=str2func(['@(t) ' num2str(double(qd))]);
            end

            if has(qdd,t)
                qdd=matlabFunction(qdd);
            else
                qdd=str2func(['@(t) ' num2str(double(qdd))]);
            end
            uqt{i_qj+5}={q;qd;qdd};
        end
        i_qj=i_qj+6;
        
    else
        error('Check joint actuation');
    end
end
end

