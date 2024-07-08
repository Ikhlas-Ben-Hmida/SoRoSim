%Function for the static equilibrium function of the linkage
%Last modified by Anup Teejo Mathew 02.03.2022
function [q,u]=statics(Tr,qu0,magnifier)

n_jact             = Tr.n_jact;
nact               = Tr.nact;
uq                 = zeros(nact,1);
if Tr.Actuated
    WrenchControlled  = Tr.WrenchControlled;

    if ~Tr.CAS
        uq(1:n_jact) = InputJointUQ0(Tr);

        for i=n_jact+1:nact
            prompt = ['Enter the strength of the soft actuator ',num2str(i-n_jact),' (N):'];
            answer = input(prompt, 's');
            uq(i)  = str2num(answer);
        end
    end
end

%initial guess definition: reference configuration
if nargin==1||isempty(qu0)
    qu0              = zeros(Tr.ndof,1);
    prompt           = "Enter the value of initial guess (normalized joint coordinates/wrenches)";
    dlgtitle         = 'Initial condition';
    definput         = {num2str(qu0')};
    opts.Interpreter = 'tex';
    answer           = inputdlg(prompt,dlgtitle,[1 75],definput,opts);
    qu0              = str2num(answer{1})';
end


disp('Time-advancing')

% qu0 = qu0./Tr.q_scale;
if nargin<=2
    magnifier = 1;
end
lsqoptions = optimoptions('lsqlin','Display','off'); %used if custom actuator strength is enabled and actuator strength is calculated using lsqlin
Func    = @(qu) Equilibrium(Tr,qu,uq,magnifier,lsqoptions);
options = optimoptions('fsolve','Algorithm','trust-region-dogleg','Display','iter','MaxFunctionEvaluations',2e10);%,'Jacobian','on'); 'trust-region-dogleg' (default), 'trust-region', and 'levenberg-marquardt'.
tic
qu      = fsolve(Func,qu0,options);
% qu = lsqnonlin(Func,qu0,[0.0001,0.001,0.001],[1,10,50]);
toc

%% plot linkage with new q
q = qu;
if Tr.Actuated
    u        = uq;
    i_jactq  = Tr.i_jactq;
    for i=1:n_jact
        if ~WrenchControlled(i)
            q(i_jactq(i)) = uq(i);
            u(i)          = qu(i_jactq(i));
        end
    end
end

% q=q.*Tr.q_scale;

%change q
save('StaticsSolution.mat','q');
if Tr.Actuated
    save('StaticsSolution.mat','q','u');
end
% end

Tr.plotq(q);
end
