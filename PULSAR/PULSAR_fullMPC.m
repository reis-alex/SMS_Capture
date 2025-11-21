addpath(genpath([pwd '\Toolbox_SMS_sym']));
clear all, clc
close all
import casadi.*

% load robot through URDFs, generate model
robot_path_pre = [pwd '\robot_isparo.urdf'];
[satelite,~] = urdf2robot_flex_visu(robot_path_pre);
[robot_pre_sim]= robot_slx_format(satelite);
satelite = SPART_casadi(robot_path_pre);

% auxiliary skew symmetric matrix
SkewSym = @(x)[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

% for quaternion integration and conversions
quaternion = quaternion();


%% System dimensions
xsat_ini = vertcat(zeros(6,1), ...                  % null initial displacements
                    zeros(3,1), ... % first 3: q_wheels,
                    zeros(5,1), ... % last 5: q_arm
                    zeros(length(satelite.idx.velocities),1));                               % null velocities

% initialize satelite's state vector
R0s = eye(3);
q0s = quaternion.rotm_to_quat(R0s); %initial quaternion, obtained directly from R0

Hs = satelite.dynamics.H;
Cs = satelite.dynamics.C;

link_positions(:,:,1) = full(satelite.kinematics.rL(R0s,xsat_ini(1:6+satelite.robot.n_q)));
pos_t = full(link_positions(:,end,1));

%% Create NDI: reaction wheels
import casadi.*

Kp1 = 1.5*eye(3);
Kv1 = Kp1*0.1;
K1  = horzcat(Kv1,Kp1);

H_sym = satelite.dynamics.H(satelite.R0,satelite.state_vars);
C_sym = satelite.dynamics.C(satelite.R0,satelite.state_vars);
virtual_control = casadi.SX.sym('v',6,1);

n0 = 6;
nq = satelite.robot.n_q;
nr = 3; 

H00 = H_sym(1:n0, 1:n0);
H0q = H_sym(1:n0, n0+1 : n0+nq);
Hqq = H_sym(n0+1 : n0+nq, n0+1 : n0+nq);

vel_sym = satelite.state_vars(satelite.idx.velocities); 
C0 = C_sym(1:n0, :) * vel_sym;            
Cq = C_sym(n0+1 : n0+nq, :) * vel_sym;     

eps_reg = 1e-9;
Hqq_reg_inv = inv(Hqq + eps_reg*eye(nq));

G = mtimes(H0q, Hqq_reg_inv);           
S = H00 - mtimes(H0q, mtimes(Hqq_reg_inv, H0q'));   

r = - C0 + mtimes(G, Cq) - mtimes(S, virtual_control); % does not depend on r0, only R0, q and qdot

% NDI for wheels
G_r = G(:, 1:nr);
lambda = 1e-8;
tau_r = mtimes((G_r' * G_r + lambda*eye(nr)) \ (G_r' ), r);

NDI_wheels = Function('NDI_fun', {satelite.R0, satelite.state_vars, virtual_control}, {tau_r}, {'R0','x','v_des'}, {'tau_r'});
NDI_compensation = casadi.Function('NDI_comp',{satelite.R0, satelite.state_vars},{mtimes((G_r' * G_r + lambda*eye(nr)) \ (G_r' ),- C0 + mtimes(G, Cq))});
NDI_PD = casadi.Function('NDI_comp',{satelite.R0, satelite.state_vars, virtual_control},{mtimes((G_r' * G_r + lambda*eye(nr)) \ (G_r' ),- mtimes(S, virtual_control))});

%% MPC
tau = casadi.SX.sym('torque_arm',satelite.robot.n_q,1);
% Define MPC problem
opt.model.states   = [satelite.state_vars];
f_mpc = satelite.dynamics.ddX(eye(3),satelite.state_vars,tau);
f_mpc = casadi.Function('fmpc',{satelite.R0,satelite.state_vars,tau},{f_mpc});

% texp = [satelite.state_vars(satelite.idx.velocities);f_mpc(satelite.R0,satelite.state_vars,tau)];
% texp = casadi.substitute(texp, satelite.state_vars(15:20),  -H00\H0q*satelite.state_vars(21:end));
% texp = casadi.substitute(texp, satelite.R0,  eye(3));
% fexp = casadi.Function('f',{satelite.state_vars([1:14 21:28]),tau},{texp([1:14 21:end])});

opt.model.function = casadi.Function('fmpc',{satelite.state_vars,tau},{vertcat(satelite.state_vars(satelite.idx.velocities), ... 
                                             f_mpc(eye(3),satelite.state_vars,tau))});

opt.model.controls = tau;
opt.continuous_model.integration = 'euler';

opt.dt          = 0.1;
opt.n_controls  = satelite.robot.n_q;          
opt.n_states    = length(opt.model.states);
opt.N           = 5;


R       = eye(satelite.robot.n_q);

opt.parameters.name{1} = 'target';
opt.parameters.name{2} = 'itm_target';
opt.parameters.dim = vertcat([3,1], [3,1]);
% opt.parameters.dim = vertcat([3,1]);

% Cost stage function
% ee_fun = satelite.kinematics.rL(casadi.SX.eye(3),satelite.state_vars);
ee_fun = satelite.kinematics.rL(casadi.SX.eye(3),vertcat(satelite.state_vars(1:14,1)));
ee_fun = casadi.Function('ee',{opt.model.states(1:6+satelite.robot.n_q)},{ee_fun(:,end)});

opt.costs.stage.parameters = opt.parameters.name(2);
opt.costs.stage.function   = @(x,u,varargin) 10*sum((ee_fun(x(1:6+satelite.robot.n_q))-varargin{:}).^2) + u'*R*u + x(1:3)'*1e5*x(1:3);% + x(15:17)'*1*x(15:17);

opt.costs.general.parameters = opt.parameters.name;
opt.costs.general.function   = @(x,u,varargin) 1000*(varargin{end}-varargin{end-1})'*(varargin{end}-varargin{end-1});

opt.constraints.states.upper  = vertcat( inf*ones(3,1),  inf*ones(3,1),  inf*ones(satelite.robot.n_q,1),  inf*pi*ones(3,1),  inf*ones(3,1),  inf*ones(3,1),  0.0939*ones(satelite.robot.n_q-3,1));
opt.constraints.states.lower  = vertcat(-inf*ones(3,1), -inf*ones(3,1), -inf*ones(satelite.robot.n_q,1), -inf*pi*ones(3,1), -inf*ones(3,1), -inf*ones(3,1), -0.0939*ones(satelite.robot.n_q-3,1));

opt.constraints.control.upper = vertcat(0.175*ones(3,1),50*ones(5,1));
opt.constraints.control.lower = -opt.constraints.control.upper;

opt.constraints.general.parameters  = opt.parameters.name(2);
opt.constraints.general.function{1} = @(x,varargin) ee_fun(x(1:6+satelite.robot.n_q))-varargin{1};
opt.constraints.general.elements{1} = 'end';
opt.constraints.general.type{1} = 'equality';

opt.constraints.general.parameters  = opt.parameters.name(2);
opt.constraints.general.function{2} = @(x,varargin) x(1:3);
opt.constraints.general.elements{2} = 'end';
opt.constraints.general.type{2} = 'equality';


opt.constraints.parameters.name  = opt.parameters.name(2);
opt.constraints.parameters.upper = inf*ones(3,1);
opt.constraints.parameters.lower = -inf*ones(3,1);

% Define inputs to optimization
opt.input.vector = opt.parameters.name(1);
opt.solver = 'ipopt';
[solver_mpc,args_mpc] = build_mpc(opt);
ine_wheels = {satelite.robot.links.inertia};
ine_wheels = ine_wheels{2};
%% simulation loop
% simulation parameters
mpc_x0   = zeros(1,length(args_mpc.vars{1}));

dq0 = casadi.Function('dq0', {satelite.R0,satelite.state_vars([1:14 21:end])}, {-H00\H0q*satelite.state_vars(21:end)});
xsat = xsat_ini;
for k = 1:20000
    k
    link_positions(:,:,k) = full(satelite.kinematics.rL(R0s(:,:,k),vertcat(xsat(1:6+satelite.robot.n_q,k))));

    mpc_input = vertcat(xsat(:,k), link_positions(:,end,1)-[1;1;0]);
    tic
    sol = solver_mpc('x0', mpc_x0, 'lbx', args_mpc.lbx, 'ubx', args_mpc.ubx, ...
                     'lbg', args_mpc.lbg, 'ubg', args_mpc.ubg, 'p', mpc_input);
    tsol(k) = toc;
    mpc_x0        = full(sol.x);
    torque(:,k) = full(sol.x(args_mpc.vars{3}));
    xsat(:,k+1) = xsat(:,k) + full(opt.dt*opt.model.function(xsat(:,k),torque(:,k))); %vertcat(xsat(satelite.idx.velocities,k),full(feval));

    % update satelite quaternion, R0 according to next omega
    [q0s(:,k+1), R0s(:,:,k+1)]  = quaternion.integrate(xsat(satelite.idx.omega0,k),q0s(:,k),opt.dt);

    itmt(:,k) = full(sol.x(end-2:end));

    % test
    tt(:,k) = full(dq0(R0s(:,:,k),xsat([1:14 21:end],k)) - xsat(15:20,k));

end

%%
figure
plot(itmt(1,:),'--r')
hold on
plot(squeeze(link_positions(1,end,:)),'-r')
plot(itmt(2,:),'--g')
plot(squeeze(link_positions(2,end,:)),'-g')
plot(itmt(3,:),'--b')
plot(squeeze(link_positions(3,end,:)),'-b')
xlabel('Time [s]')
ylabel('Position [m]')
grid on
title('Position of end effector')

figure
plot(xsat(1,:),'r')
hold on
plot(xsat(2,:),'g')
plot(xsat(3,:),'b')
xlabel('Time [s]')
ylabel('Position [rad]')
grid on
title('Angular position of the base')


figure
plot(xsat(21,:),'r')
hold on
plot(xsat(22,:),'g')
plot(xsat(23,:),'b')
xlabel('Time [s]')
ylabel('Velocity [rad/s]')
grid on
title('Angular velocity of the reaction wheels')

figure
plot(torque(1,:),'r')
hold on
plot(torque(2,:),'g')
plot(torque(3,:),'b')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on
title('Reaction wheels torques')

figure
plot(torque(4,:),'r')
hold on
plot(torque(5,:),'g')
plot(torque(6,:),'b')
plot(torque(7,:),'c')
plot(torque(8,:),'m')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on
title('Manipulator joint torques')