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
                    20*pi/180,30*pi/180,20*pi/180,10*pi/180,25*pi/180, ... % last 5: q_arm
                    zeros(length(satelite.idx.velocities),1));                               % null velocities

% used as reference
% xsat_ini = vertcat(zeros(6,1), ...                  % null initial displacements
%                     zeros(3,1), ... % first 3: q_wheels,
%                     225*pi/180,150*pi/180,-50*pi/180,30*pi/180,25*pi/180, ... % last 5: q_arm
%                     zeros(length(satelite.idx.velocities),1)); 

% initialize satelite's state vector
R0s = eye(3);
q0s = quaternion.rotm_to_quat(R0s); %initial quaternion, obtained directly from R0

Hs = satelite.dynamics.H;
Cs = satelite.dynamics.C;

link_positions(:,:,1) = full(satelite.kinematics.rL(R0s,xsat_ini(1:6+satelite.robot.n_q)));
pos_t = full(link_positions(:,end,1));

% since r0 = 0, the reference computed here is w.r.t. the body frame

plot3(link_positions(1,7:end,1),link_positions(2,7:end,1),link_positions(3,7:end,1),'-og')
hold on
plot3(link_positions(1,end,1),link_positions(2,end,1),link_positions(3,end,1),'-xr','LineWidth',10)
plot3(link_positions(1,7,1),link_positions(2,7,1),link_positions(3,7,1),'or','LineWidth',10)
axis auto

%%
H_sym = satelite.dynamics.H(satelite.R0,satelite.state_vars);
C_sym = satelite.dynamics.C(satelite.R0,satelite.state_vars);

n0 = 6;
nq = satelite.robot.n_q;
nr = 3; 

H00 = H_sym(1:n0, 1:n0);
H0q = H_sym(1:n0, n0+1 : n0+nq);
Hqq = H_sym(n0+1 : n0+nq, n0+1 : n0+nq);

moment_satelite = casadi.Function('h_sat',{satelite.R0,satelite.state_vars},{H00*satelite.state_vars(satelite.idx.velocities(1:6))});
moment_wheels = casadi.Function('h_wheel',{satelite.R0,satelite.state_vars},{H0q(:,1:3)*satelite.state_vars(satelite.idx.velocities(7:9))});
moment_arm = casadi.Function('h_arm',{satelite.R0,satelite.state_vars},{H0q(:,4:end)*satelite.state_vars(satelite.idx.velocities(10:14))});

% R0 function
angles = satelite.state_vars(1:3);
yaw = angles(3);
pitch = angles(2);
roll = angles(1);
c1 = cos(yaw);
s1 = sin(yaw);
c2 = cos(pitch);
s2 =  sin(pitch);
c3 = cos(roll);
s3 =  sin(roll);
R = [c1*c2, c1*s2*s3 - s1*c3, c1*s2*c3 + s1*s3;
     s1*c2, s1*s2*s3 + c1*c3, s1*s2*c3 - c1*s3;
     -s2,   c2*s3,             c2*c3];
Rf = casadi.Function('R0',{satelite.state_vars(1:3)},{R});

%% Create NDI: reaction wheels
import casadi.*
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

r = - C0 + mtimes(G, Cq) - mtimes(S, virtual_control);

% NDI for wheels
G_r = G(:, 1:nr);
lambda = 1e-8;
tau_r = mtimes((G_r' * G_r + lambda*eye(nr)) \ (G_r' ), r);

Kp1 = 50*eye(3);
Kv1 = 0.01*eye(3);
K1  = horzcat(Kp1,Kv1);

theta_d = casadi.SX.sym('theta_d',3,1);
dtheta_d = casadi.SX.sym('dtheta_d',3,1);

v = vertcat(K1*vertcat((satelite.state_vars(1:3)-theta_d),satelite.state_vars(15:17)-dtheta_d), zeros(3,1));

tau_r = casadi.substitute(tau_r, virtual_control,v);
tau_r = casadi.substitute(tau_r, satelite.state_vars(15:20),  -H00\H0q*satelite.state_vars(21:end));
tau_r = casadi.substitute(tau_r, satelite.R0,  R);

NDI_wheels = casadi.Function('NDI_fun', {satelite.state_vars([1:14 21:end]),theta_d,dtheta_d}, {tau_r}, {'x','target','dtarget'}, {'tau_r'});

%% MPC
tau = casadi.SX.sym('torque_arm',satelite.robot.n_q,1);
desired0 = casadi.SX.sym('desired0',6,1);

% Define MPC problem
opt.model.states   = [satelite.state_vars([1:14 21:end])];
f_mpc = vertcat(satelite.dynamics.ddX(satelite.R0,satelite.state_vars,tau));
f_mpc = casadi.substitute(f_mpc, satelite.state_vars(15:20),  -H00\H0q*satelite.state_vars(21:end));

texp = [satelite.state_vars(satelite.idx.velocities);f_mpc(7:end);];
texp = casadi.substitute(texp, satelite.state_vars(15:20),  -H00\H0q*satelite.state_vars(21:end));
texp = casadi.substitute(texp, satelite.R0,  R);
texp = casadi.substitute(texp,tau(1:3),NDI_wheels(satelite.state_vars([1:14 21:end]),desired0(1:3),desired0(4:6)));

q0dot = -H00\H0q*satelite.state_vars(21:end);
omega0 = q0dot(1:3);
omega0 = casadi.substitute(omega0, satelite.R0,  R);
omega0 = casadi.Function('omega0',{satelite.state_vars([1:14 21:end])},{omega0});

accel =  casadi.substitute(f_mpc(1:3),satelite.state_vars(15:20),  -H00\H0q*satelite.state_vars(21:end));
accel = casadi.substitute(accel, satelite.R0,  R);
accel = casadi.substitute(accel,tau(1:3),NDI_wheels(satelite.state_vars([1:14 21:end]),desired0(1:3),desired0(4:6)));
accel = casadi.Function('accel',{satelite.state_vars([1:14 21:end]), tau(4:end), desired0(1:3),desired0(4:6)},{accel});

opt.model.function = casadi.Function('fmpc',{vertcat(satelite.state_vars([1:14 21:28])),vertcat(tau(4:end),desired0(1:3),desired0(4:6))},{texp});

% opt.model.controls = tau(4:end);
opt.model.controls = vertcat(tau(4:end),desired0);
opt.continuous_model.integration = 'euler';

opt.dt          = 0.1;
opt.n_controls  = length(opt.model.controls);          
opt.n_states    = length(opt.model.states);
opt.N           = 4;

% R       = blkdiag(1*eye(3),1*eye(5));
p = 10*pi/180;
P = eye(3)*1/(p^2);

opt.parameters.name{1} = 'target';
opt.parameters.name{2} = 'itm_target';
opt.parameters.name{3} = 'aux';
opt.parameters.dim = vertcat([3,1], [3,1], [1,1]);

% Cost stage function
ee_fun = satelite.kinematics.rL(R,vertcat(satelite.state_vars(1:14,1)));
ee_fun = casadi.Function('ee',{opt.model.states(1:6+satelite.robot.n_q)},{R'*ee_fun(:,end)});

opt.costs.stage.parameters = opt.parameters.name(1:3);
opt.costs.stage.function   = @(x,u,varargin) sum((ee_fun(x(1:6+satelite.robot.n_q))-varargin{:}(4:6)).^2) + varargin{:}(1)*(x(1:3,end)'*1e6*x(1:3,end) + (x(15:17,end)-10)'*1e2*(x(15:17,end)-10));% + varargin{:}(7)*(omega0(x)'*0*omega0(x) + 1e5*x(1:3)'*x(1:3));% + x(18:22)'*x(18:22))*1e3;%u(6:8)'*u(6:8)*1e2 +

opt.costs.general.parameters = opt.parameters.name(1:2);
opt.costs.general.function   = @(x,u,varargin) 1e5*(varargin{2}-varargin{1})'*(varargin{2}-varargin{1});

opt.constraints.states.upper  = vertcat( inf*ones(3,1),  inf*ones(3,1), inf*ones(3,1),  inf*ones(satelite.robot.n_q-3,1),  3600*2*pi/60*ones(3,1),  0.09*ones(satelite.robot.n_q-3,1));
opt.constraints.states.lower  = vertcat(-inf*ones(3,1), -inf*ones(3,1), -inf*ones(3,1), -inf*ones(satelite.robot.n_q-3,1), -3600*2*pi/60*ones(3,1), -0.09*ones(satelite.robot.n_q-3,1));

opt.constraints.control.upper = vertcat(50*ones(5,1), p*ones(3,1), 0.2*ones(3,1)); %0.175*ones(3,1),
opt.constraints.control.lower = -opt.constraints.control.upper;

opt.constraints.general.parameters  = opt.parameters.name(2);
opt.constraints.general.function{1} = @(x,varargin) (ee_fun(x(1:6+satelite.robot.n_q))-varargin{:}(1:3));
opt.constraints.general.elements{1} = 'end';
opt.constraints.general.type{1}     = 'equality';

opt.constraints.general.function{2} = @(x,u,varargin) vertcat(NDI_wheels(x,u(6:8),u(9:11))-0.175*ones(3,1), ...
                                                             -NDI_wheels(x,u(6:8),u(9:11))-0.175*ones(3,1));
opt.constraints.general.elements{2} = 'N';
opt.constraints.general.type{2} = 'inequality';


opt.constraints.parameters.name  = opt.parameters.name(2);
opt.constraints.parameters.upper =  vertcat(inf*ones(3,1));
opt.constraints.parameters.lower = -vertcat(inf*ones(3,1));

% Define inputs to optimization
opt.input.vector = opt.parameters.name([1 3]);
opt.solver = 'ipopt';
[solver_mpc,args_mpc] = build_mpc(opt);
ine_wheels = {satelite.robot.links.inertia};
ine_wheels = ine_wheels{2};


%% simulation loop
% simulation parameters
mom_arm = 0.01*ones(6,1);
mom_wheels = 0.01*ones(6,2);
mpc_x0   = zeros(1,length(args_mpc.vars{1}));

dq0 = casadi.Function('dq0', {satelite.R0,satelite.state_vars([1:14 21:end])}, {-H00\H0q*satelite.state_vars(21:end)});
xsat(:,1) = vertcat(xsat_ini([1:14 21:end],1));
sw = 0;
sw2 = 1;
feval = zeros(31,1);
pred_states = zeros(opt.n_states,opt.N);
%%
for k = 1:2000
    k

    link_positions(:,:,k) = full(R0s(:,:,k)'*satelite.kinematics.rL(R0s(:,:,k),vertcat(xsat(1:6+satelite.robot.n_q,k))));
%     ref =  [2.806; -0.9745; -0.1926];
ref =  [5.0528; 0.6535; -2.0646];

    if k<50 %norm(link_positions(:,end,k)-ref)>=0.0005 && sw == 0
        htg = zeros(3,1);
        auxi = 0;
        sig = zeros(3,1);
    else
        if sw2 == 1
            auxi = 1;
            sw = 1;
            sw2 = 0;
        end
    end 

    mpc_input = vertcat(xsat(:,k),ref,auxi);

    sol = solver_mpc('x0', mpc_x0, 'lbx', args_mpc.lbx, 'ubx', args_mpc.ubx, ...
        'lbg', args_mpc.lbg, 'ubg', args_mpc.ubg, 'p', mpc_input);
    mpc_x0        = full(sol.x);

    % gather torques
    torque_arm(:,k) = full(sol.x(args_mpc.vars{3}(1:5)));
    desired(:,k) = full(sol.x(args_mpc.vars{3}(6:end)));
    torque_wheels(:,k) = full(NDI_wheels(xsat(:,k),desired(1:3,k),desired(4:6,k)));

    pred_states = reshape(full(sol.x(1:length(opt.model.states)*opt.N)),length(opt.model.states),opt.N);
    feval = full(opt.model.function(xsat(:,k),vertcat(torque_arm(:,k),desired(:,k))));
    xsat(:,k+1) = xsat(:,k) + opt.dt*feval;

    check_feas(solver_mpc.stats())

    % update satelite quaternion, R0 according to next omega
    tt(:,k) = full(dq0(R0s(:,:,k),xsat(1:22,k)));
    [q0s(:,k+1), R0s(:,:,k+1)]  = quaternion.integrate(tt(1:3,k),q0s(:,k),opt.dt);
    xsat_full(:,k) = vertcat(xsat(1:14,k),tt(:,k),xsat(15:22,k));

    % compute instantaneous momenta
    mom_arm(:,k) = full(moment_arm(R0s(:,:,k),xsat_full(:,k)));
    mom_sate(:,k) = full(moment_satelite(R0s(:,:,k),xsat_full(:,k)));
    mom_wheels(:,k) = full(moment_wheels(R0s(:,:,k),xsat_full(:,k)));

    itmt(:,k) = full(sol.x(end-2:end));
    accelerations(:,k) = full(accel(xsat(:,k),torque_arm(:,k),desired(1:3,k),desired(4:6,k)));
    cost(k) = full(sol.f);
end

%%
close all


figure
plot(itmt(1,:),'-.r')
hold on
plot(ref(1,:)*ones(1,length(squeeze(link_positions(1,end,:)))),'--r')
plot(squeeze(link_positions(1,end,:)),'-r')
plot(itmt(2,:),'-.g')
plot(ref(2,:)*ones(1,length(squeeze(link_positions(1,end,:)))),'--g')
plot(squeeze(link_positions(2,end,:)),'-g')
plot(itmt(3,:),'-.b')
plot(ref(3,:)*ones(1,length(squeeze(link_positions(1,end,:)))),'--b')
plot(squeeze(link_positions(3,end,:)),'-b')
xlabel('Time [s]')
ylabel('Position [m]')
grid on
title('Position of end effector')
legend('EE_x','EE_x^{ref}','EE_y','EE_y^{ref}','EE_z','EE_z^{ref}','interpreter','latex')

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
plot(xsat(15,:),'r')
hold on
plot(xsat(16,:),'g')
plot(xsat(17,:),'b')
xlabel('Time [s]')
ylabel('Velocity [rad/s]')
grid on
title('Angular velocity of the reaction wheels')

figure
plot(xsat(18,:),'r')
hold on
plot(xsat(19,:),'g')
plot(xsat(20,:),'b')
plot(xsat(21,:),'c')
plot(xsat(22,:),'m')
xlabel('Time [s]')
ylabel('Velocity [rad/s]')
grid on
title('Angular velocity of the manipulator joints')

figure
plot(torque_wheels(1,:));
hold on;
plot(torque_wheels(2,:))
plot(torque_wheels(3,:))
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on
title('Reaction wheels torques')

figure
stairs(torque_arm(1,:),'r')
hold on
stairs(torque_arm(2,:),'g')
stairs(torque_arm(3,:),'b')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on
title('Manipulator joint torques')

figure
subplot(131)
plot(mom_sate(1,:),'b'); hold on; plot(mom_arm(1,:),'-g'),plot(mom_wheels(1,:),'--r')
subplot(132)
plot(mom_sate(2,:),'b'); hold on; plot(mom_arm(2,:),'-g'),plot(mom_wheels(2,:),'--r')
subplot(133)
plot(mom_sate(3,:),'b'); hold on; plot(mom_arm(3,:),'-g'),plot(mom_wheels(3,:),'--r')
legend('mom sat','mom bras', 'mom roues')

figure
plot(mom_wheels(1,:),'-r')
hold on
plot(mom_wheels(2,:),'-g')
plot(mom_wheels(3,:),'-b')

figure
subplot(121)
plot(desired(1,:),'r'); 
hold on;
plot(desired(2,:),'g')
plot(desired(3,:),'b')
subplot(122)
plot(desired(4,:),'r'); 
hold on;
plot(desired(5,:),'g')
plot(desired(6,:),'b')
%%
figure
hold on
plot3(link_positions(1,7,1),link_positions(2,7,1),link_positions(3,7,1),'or','LineWidth',10)
hold on
for i = 1:k
    plot3(link_positions(1,7:end,i),link_positions(2,7:end,i),link_positions(3,7:end,i),'-og')
    plot3(link_positions(1,end,i),link_positions(2,end,i),link_positions(3,end,i),'-xr','LineWidth',10)
    axis auto
    view([-37.5 30])
    drawnow
    grid
end

%%
function [] = check_feas(stats)
        if ~stats.success
            disp('Solver failed.');
            disp(['Return status: ', stats.return_status]);

            if contains(stats.return_status, 'Infeasible')
                error('Problem is infeasible!');
            end
        end
end