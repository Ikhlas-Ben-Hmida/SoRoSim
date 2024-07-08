%Function to calculate actuator strength as a function of system dynamics
%Last modified by Anup Teejo Mathew 02.03.2022

function u = CustomActuatorStrength(Tr,q,g,J,t,qd,eta,Jd,M,C,F,Bq)


%Tr: Linkage element,
%q and qd: joint coordinates and their time derivatives,
%g, J, Jd, and eta: transformation matrix, Jacobian, time derivative of jacobian, and screw velocity at every significant point of the linkage
%t:  time
%M,C,F,Bq: generalized mass, coriolis, force, and actuation matrices.
%u should be (nactx1) column vector where nact is the total number of actuators.


% global us ts es
control_type   = 6;   % 3 for 3D, or 6 for 6D control
control_target = 'T'; % 'P' for desired point, or 'T' for desired trajectory control

%% Desired Tip Trajectory / position and orientation

if control_target == 'T' % Controller follows a set trajectory

    w = pi/8;
    
    g_fk = [0.8094   -0.5872    0    0.9574;...
            0.5872    0.8094    0    0.2162;...
            0         0         1         0;...
            0         0         0    1.0000];
    
    R               = [1 0 0;0 cos(w*t) -sin(w*t);0 sin(w*t) cos(w*t)];
    Rdot            = [0 0 0;0 -w*sin(w*t) -w*cos(w*t);0 w*cos(w*t) -w*sin(w*t)];
    Rdotdot         = [0 0 0;0 -w^2*cos(w*t) w^2*sin(w*t);0 -w^2*sin(w*t) -w^2*cos(w*t)];
    g_bar           = [R,zeros(3,1);0,0,0,1]*g_fk*[R',zeros(3,1);0,0,0,1];
    g_bar_dot       = [Rdot,zeros(3,1);0,0,0,0]*g_fk*[R',zeros(3,1);0,0,0,1]+[R,zeros(3,1);0,0,0,1]*g_fk*[Rdot',zeros(3,1);0,0,0,0];
    g_bar_dotdot    = [Rdotdot,zeros(3,1);0,0,0,0]*g_fk*[R',zeros(3,1);0,0,0,1]+2*[Rdot,zeros(3,1);0,0,0,0]*g_fk*[Rdot',zeros(3,1);0,0,0,0]+[R,zeros(3,1);0,0,0,1]*g_fk*[Rdotdot',zeros(3,1);0,0,0,0];
    

elseif control_target == 'P' % Controller finds a desired point

     
g_bar = [-0.6017   -0.2149    0.7693    0.1468
    0.7349    0.2283    0.6386    0.5236
   -0.3129    0.9496    0.0206   -0.6633
         0         0         0    1.0000];
    g_bar_dot    = zeros(4,4);
    g_bar_dotdot = zeros(4,4);
    
end

%% Tip position and orientation
g_tip   = g(end-3:end,:);
g_error = ginv(g_tip)*g_bar;

%% Tip velocity and acceleration
eta_tip      = eta(end-5:end,:);

%% Desired Tip Velocity and Accleration (in the target frame)
eta_bar_hat     = (g_bar^-1)*g_bar_dot;
eta_bar_hat_dot = -eta_bar_hat^2+g_bar^-1*g_bar_dotdot;

eta_bar     = [eta_bar_hat(3,2);eta_bar_hat(1,3);eta_bar_hat(2,1);eta_bar_hat(1,4);eta_bar_hat(2,4);eta_bar_hat(3,4)];
eta_bar_dot = [eta_bar_hat_dot(3,2);eta_bar_hat_dot(1,3);eta_bar_hat_dot(2,1);eta_bar_hat_dot(1,4);eta_bar_hat_dot(2,4);eta_bar_hat_dot(3,4)];

%% Error in tip velocity and Accelaration (in tip frame)
e_eta       = dinamico_Adjoint(g_error)*eta_bar-eta_tip;                     %in tip frame
eta_bar_dot = dinamico_adj(e_eta)*dinamico_Adjoint(g_error)*eta_bar+... %in tip frame
              dinamico_Adjoint(g_error)*eta_bar_dot;

%% Tip jacobian and it's derivative
J_tip       = J(end-5:end,:);
Jd_tip      = Jd(end-5:end,:);

%% Dynamic Equation Coeffecients
D       = Tr.D;
K       = Tr.K;

%% Controller
if control_type == 3 %3D
    %%
    J_tip_pos       = J_tip(4:6,:);
    Jd_tip_pos      = Jd_tip(4:6,:);
    e_pos           = g_error(1:3,4);
    e_eta_pos       = e_eta(4:6);
    eta_bar_dot_pos = eta_bar_dot(4:6);

    m_n         = size(J_tip_pos);
    full_rank   = min(m_n);
    m           = rank(J_tip_pos);

    if m ~= full_rank
        fprintf('Rank = %d \n',m)
        fprintf('\n Jacobian is not full rank\n')
%         pause
    end
    K_p = 120*eye(3);
    K_d = 120*eye(3);
    
    options = optimoptions('lsqlin','Algorithm','interior-point','Display','off'); % slows down, fix it by passing as argument
    u = lsqlin(J_tip_pos*(M^-1)*Bq,eta_bar_dot_pos+K_d*(e_eta_pos)+(K_p*(e_pos))+J_tip_pos*(M^-1)*((C+D)*qd+K*q-F)-Jd_tip_pos*qd,[],[],[],[],-200*ones(Tr.n_sact,1),0*ones(Tr.n_sact,1),[],options);


elseif control_type == 6 %6D

    e = piecewise_logmap(g_error);

    K_p        = 100*eye(6);
    K_d        = 100*eye(6);

    lsqoptions = optimoptions('lsqlin','Display','off'); %is defined in dynamics.m
    u = lsqlin(J_tip*(M^-1)*Bq,eta_bar_dot+K_d*(e_eta)+(K_p*(e))+J_tip*(M^-1)*((C+D)*qd+K*q-F)-Jd_tip*qd,[],[],[],[],-50*ones(Tr.n_sact,1),0.001*ones(Tr.n_sact,1),[],lsqoptions);
%     u = lsqnonneg(J_tip*(M^-1)*Bq,eta_bar_dot+K_d*(e_eta)+(K_p*(e))+J_tip*(M^-1)*((C+D)*qd+K*q-F)-Jd_tip*qd,[],lsqoptions);
e
% us=[us;u'];
% ts=[ts;t];
% es=[es;e'*e];

end

