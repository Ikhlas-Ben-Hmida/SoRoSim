%Function for the dynamic simulation of the linkage
%Last modified by Anup Teejo Mathew 23/05/2021
function [t,qqd] = dynamics(S)

nact  = S.nact;
uqt   = cell(nact,1);
ndof  = S.ndof;

if S.Actuated
    
    n_jact             = S.n_jact;
    i_jactq            = S.i_jactq;
    WrenchControlled   = S.WrenchControlled;
    
    if ~S.CAS
        
        uqt(1:n_jact) = InputJointUQt(S);

        for i=n_jact+1:nact

                prompt = ['Enter the strength (N) of the soft actuator ',num2str(i-S.n_jact),' as a function of t (time)'...
                          '\n[Examples: -10-5*t, -50*t+(50*t-50)*heaviside(t-1), 50*sin(2*pi*t)]: '];
                funstr = input(prompt, 's');
                uqt{i} = str2func( ['@(t) ' funstr ] );

        end
    end
    
end


%initial guess definition: reference configuration
q0                 = zeros(1,ndof);
ndof               = S.ndof;
qd0                = zeros(1,ndof);

if S.Actuated
    i_qControlled      = i_jactq.*~WrenchControlled;
    i_qControlled      = i_qControlled(i_qControlled~=0);
    q0(i_qControlled)  = [];
    qd0(i_qControlled) = [];
end

prompt           = {'Enter the value of q_{0}','qdot_{0}','Simulation time (s)'};
dlgtitle         = 'Initial condition';
definput         = {num2str(q0),num2str(qd0),'5'};
opts.Interpreter = 'tex';
answer           = inputdlg(prompt,dlgtitle,[1 75],definput,opts);

q0   = str2num(answer{1})';
qd0  = str2num(answer{2})';
tmax = str2num(answer{3});

if S.Actuated
    q0_pass          = zeros(ndof,1);
    qd0_pass         = zeros(ndof,1);
    i                = 1:ndof;
    i(i_qControlled) = [];
    q0_pass(i)       = q0;
    qd0_pass(i)      = qd0;
    q0               = q0_pass;
    qd0              = qd0_pass;

    for i=1:n_jact
        if ~WrenchControlled(i)
            q0(i_jactq(i))   = uqt{i}{1}(0);
            qd0(i_jactq(i))  = uqt{i}{2}(0);
        end
    end
end

qqd0 = [q0./S.q_scale ; qd0./S.q_scale];

profile on
tic
[t,qqd] = ode45(@(t,qqd) S.derivatives(t,qqd,uqt),[0 tmax],qqd0);
toc
profile off

qqd = qqd.*repmat(S.q_scale',length(t),2);

save('DynamicsSolution.mat','t','qqd');

quest  = 'Generate output video of the simulation?';
Answer = questdlg(quest,'Plot Output','Yes','No','Yes');

switch Answer
    case 'Yes'
        plotqqd(S,t,qqd);
end

end

