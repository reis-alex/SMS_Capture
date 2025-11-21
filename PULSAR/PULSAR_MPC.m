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
xsat = vertcat(zeros(6,1), ...                  % null initial displacements
    zeros(3,1), ... % first 3: q_wheels,
    zeros(5,1), ... % last 5: q_arm
    zeros(length(satelite.idx.velocities),1));                               % null velocities

% initialize satelite's state vector
R0s = eye(3);
q0s = quaternion.rotm_to_quat(R0s); %initial quaternion, obtained directly from R0

Hs = satelite.dynamics.H;
Cs = satelite.dynamics.C;

link_positions(:,:,1) = full(satelite.kinematics.rL(R0s,xsat));
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
tau_m = casadi.SX.sym('torque_arm',satelite.robot.n_q-3,1);
% Define MPC problem
opt.model.states   = [satelite.state_vars];
f_mpc = satelite.dynamics.ddX(eye(3),vertcat(zeros(6,1),satelite.state_vars(7:end)),vertcat(zeros(3,1),tau_m));

% omega0 is the NDI, the rest remains the same
f_mpc = vertcat(-K1*vertcat(satelite.state_vars(satelite.idx.positions(1:3)), ...
                          satelite.state_vars(satelite.idx.omega0)), ...
                f_mpc(4:end));

opt.model.function = [satelite.state_vars(satelite.idx.velocities);f_mpc];
opt.model.controls = tau_m;
opt.continuous_model.integration = 'euler';

opt.dt          = 0.1;
opt.n_controls  = satelite.robot.n_q-3;          
opt.n_states    = length(opt.model.states);
opt.N           = 5;


R       = eye(satelite.robot.n_q-3);

opt.parameters.name{1} = 'target';
% opt.parameters.name{2} = 'itm_target';
% opt.parameters.dim = vertcat([3,1], [3,1]);
opt.parameters.dim = vertcat([3,1]);

% Cost stage function
ee_fun = satelite.kinematics.rL(casadi.SX.eye(3),satelite.state_vars);
ee_fun = casadi.Function('ee',{satelite.state_vars},{ee_fun(:,end)});

opt.costs.stage.parameters = opt.parameters.name(1);
opt.costs.stage.function   = @(x,u,varargin) 10*sum((ee_fun(x)-varargin{:}(1:3)).^2) + u'*R*u;% + 100*(mu_b0(eye(3),zeros(3,1),x(7:11),zeros(3,1))-0.15)^2;% - 100*(mu_b0(x(7:12)))^2;
% opt.costs.stage.function   = @(x,u,varargin) 10*sum((x(10:14)-varargin{:}).^2) + u'*R*u;

% opt.costs.general.parameters = opt.parameters.name;
% opt.costs.general.function   = @(x,u,varargin) 1000*(varargin{end}-varargin{end-1})'*(varargin{end}-varargin{end-1});

opt.constraints.states.upper  = vertcat( inf*ones(3,1),  inf*ones(3,1),  inf*ones(satelite.robot.n_q,1),  inf*pi*ones(3,1),  inf*ones(3,1),  inf*ones(3,1),  0.0939*ones(satelite.robot.n_q-3,1));
opt.constraints.states.lower  = vertcat(-inf*ones(3,1), -inf*ones(3,1), -inf*ones(satelite.robot.n_q,1), -inf*pi*ones(3,1), -inf*ones(3,1), -inf*ones(3,1), -0.0939*ones(satelite.robot.n_q-3,1));

opt.constraints.control.upper = vertcat(50*ones(satelite.robot.n_q-3,1));
opt.constraints.control.lower = -opt.constraints.control.upper;

opt.constraints.general.parameters  = opt.parameters.name(1);
opt.constraints.general.function{1} = @(x,varargin) ee_fun(x)-varargin{1};
opt.constraints.general.elements{1} = 'end';
opt.constraints.general.type{1} = 'equality';

% opt.constraints.general.parameters  = opt.parameters.name(1);
% opt.constraints.general.function{1} = @(x,varargin) x(10:14)-varargin{1};
% opt.constraints.general.elements{1} = 'end';
% opt.constraints.general.type{1} = 'equality';

% opt.constraints.parameters.name  = opt.parameters.name(2);
% opt.constraints.parameters.upper = inf*ones(3,1);
% opt.constraints.parameters.lower = -inf*ones(3,1);

% Define inputs to optimization
opt.input.vector = opt.parameters.name([1]);
opt.solver = 'ipopt';
[solver_mpc,args_mpc] = build_mpc(opt);
ine_wheels = {satelite.robot.links.inertia};
ine_wheels = ine_wheels{2};
%% simulation loop
% simulation parameters
dt = 0.01;
total_moment = [];
feval = zeros(length(satelite.state_vars)/2,1);
mpc_x0   = zeros(1,length(args_mpc.vars{1}));

for k = 1:20000
    k
    link_positions(:,:,k) = full(satelite.kinematics.rL(R0s(:,:,k),xsat(:,k)));

    % simulate free robot
    vr_wheels(:,k) = -vertcat(K1*vertcat(xsat(1:3,k)-[0;0;0.0005*(k-1)/1000],xsat(satelite.idx.omega0,k)-[0;0;0]),zeros(3,1));
    torque_wheels(:,k) = vertcat(full(NDI_wheels(R0s(:,:,k),xsat(:,k),vr_wheels(:,k))), zeros(5,1));

    mpc_input = vertcat(xsat(:,k), link_positions(:,end,1)-[.0;0.0;0]);
%     mpc_input = vertcat(xsat(:,k), [0;0;0.1;0;0]);
%     tic
%     sol = solver_mpc('x0', mpc_x0, 'lbx', args_mpc.lbx, 'ubx', args_mpc.ubx, ...
%                      'lbg', args_mpc.lbg, 'ubg', args_mpc.ubg, 'p', mpc_input);
%     tsol(k) = toc;
%     mpc_x0        = full(sol.x);
%     torque_arm(:,k) = vertcat(zeros(3,1),full(sol.x(args_mpc.vars{3})));
    torque_arm(:,k) = vertcat(zeros(3,1),zeros(5,1));

    torque(:,k) = torque_wheels(:,k)+torque_arm(:,k);
    feval = satelite.dynamics.ddX(R0s(:,:,k),xsat(:,k),torque(:,k));
    xsat(:,k+1) = xsat(:,k) + dt*vertcat(xsat(satelite.idx.velocities,k),full(feval));

    % update satelite quaternion, R0 according to next omega
    [q0s(:,k+1), R0s(:,:,k+1)]  = quaternion.integrate(xsat(satelite.idx.omega0,k),q0s(:,k),dt);
    momentum_wheels(:,k) = ine_wheels*xsat(21:23,k);
    torquew(:,k) = full(ine_wheels*feval(7:9));
    ndi_comp(:,k) = full(NDI_compensation(R0s(:,:,k),xsat(:,k)));
    ndi_pd(:,k) = full(NDI_PD(R0s(:,:,k),xsat(:,k),vr_wheels(:,k)));
end

