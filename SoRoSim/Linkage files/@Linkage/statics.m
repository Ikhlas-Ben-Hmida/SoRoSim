%Function for the static equilibrium function of the linkage
%Last modified by Anup Teejo Mathew 23/05/2021
function [q,u]=statics(S,qu0) 

n_jact             = S.n_jact;
nact               = S.nact;
uq                 = zeros(nact,1);
if S.Actuated
    WrenchControlled   = S.WrenchControlled;
    
    if ~S.CAS
        uq(1:n_jact)       = InputJointUQ0(S);

        for i=n_jact+1:nact
                prompt=['Enter the strength of the soft actuator ',num2str(i-n_jact),' (N):'];
                answer = input(prompt, 's');
                uq(i)= str2num(answer);
        end
    end
end

%initial guess definition: reference configuration
if nargin==1
    qu0              = zeros(S.ndof,1);
    prompt           = "Enter the value of initial guess (normalized joint coordinates/wrenches)";
    dlgtitle         = 'Initial condition';
    definput         = {num2str(qu0')};
    opts.Interpreter = 'tex';
    answer           = inputdlg(prompt,dlgtitle,[1 75],definput,opts);
    qu0              = str2num(answer{1})';
end

tic

disp('Time-advancing')

Func=@(qtau) Equilibrium(S,qtau,uq); 
options=optimoptions('fsolve','Display','iter');%,'Jacobian','on');
qu=fsolve(Func,qu0./S.q_scale,options);

toc
%% plot linkage with new q
q = qu;
if S.Actuated
    u        = uq;
    i_jactq  = S.i_jactq;
    for i=1:n_jact
        if ~WrenchControlled(i)
            q(i_jactq(i)) = uq(i);
            u(i)          = qu(i_jactq(i));
        end
    end
end

q=q.*S.q_scale;

%change q
save('StaticsSolution.mat','q');
if S.Actuated
    save('StaticsSolution.mat','q','u');
end

S.plotq(q);
end
