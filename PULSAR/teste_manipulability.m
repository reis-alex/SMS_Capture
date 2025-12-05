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

hsat = H0q(:,4:end)*satelite.state_vars(satelite.idx.velocities(10:14)) + H0q(:,1:3)*satelite.state_vars(satelite.idx.velocities(7:9));
hwheel = H0q(:,1:3)*satelite.state_vars(satelite.idx.velocities(7:9));
harm = H0q(:,4:end)*satelite.state_vars(satelite.idx.velocities(10:14));

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

EE_idx = find(strcmp({satelite.robot.links.name},'Link_6_'));
Jee = satelite.Jacob(EE_idx);
mu = casadi.substitute(Jee.Jm*Jee.Jm', satelite.R0,  R);
mu = casadi.Function('mu',{satelite.state_vars([1:14 21:end])},{-log(det(mu+1e-6*eye(6)))});

mu2 = casadi.substitute(Jee.Jm(:,1:3)-Jee.J0(:,1:3)*inv(satelite.robot.links(1).inertia)*H0q(1:3,1:3),satelite.R0,  R);
mu2 = mu2*mu2';
mu2 = casadi.Function('mu',{satelite.state_vars([1:14 21:end])},{-log(det(mu2+1e-6*eye(6)))});
%% MPC
tau = casadi.SX.sym('torque_arm',satelite.robot.n_q,1);

% Define MPC problem
opt.model.states   = [satelite.state_vars([1:14 21:end])];
f_mpc = vertcat(satelite.dynamics.ddX(satelite.R0,satelite.state_vars,tau));
f_mpc = casadi.substitute(f_mpc, satelite.state_vars(15:20),  -H00\H0q*satelite.state_vars(21:end));

texp = [satelite.state_vars(satelite.idx.velocities);f_mpc(7:end);];
texp = casadi.substitute(texp, satelite.state_vars(15:20),  -H00\H0q*satelite.state_vars(21:end));
texp = casadi.substitute(texp, satelite.R0,  R);

hsatexp = hsat;
hsatexp = casadi.substitute(hsatexp, satelite.R0,  R);
hsatf = casadi.Function('Hsat',{satelite.state_vars([1:14 21:end])},{hsatexp(1:3)});

harmexp = harm;
harmexp = casadi.substitute(harmexp, satelite.R0,  R);
harmf = casadi.Function('Hsat',{satelite.state_vars([1:14 21:end])},{harmexp(1:3)});

opt.model.function = casadi.Function('fmpc',{vertcat(satelite.state_vars([1:14 21:28])),tau},{texp});

opt.model.controls = tau;
opt.continuous_model.integration = 'euler';

opt.dt          = 0.1;
opt.n_controls  = satelite.robot.n_q;          
opt.n_states    = length(opt.model.states);
opt.N           = 4;

% R       = blkdiag(1*eye(3),1*eye(5));

opt.parameters.name{1} = 'target';
opt.parameters.name{2} = 'itm_target';
opt.parameters.name{3} = 'aux';
opt.parameters.name{4} = 'sign';
opt.parameters.dim = vertcat([3,1], [3,1], [1,1], [3,1]);

% Cost stage function
ee_fun = satelite.kinematics.rL(R,vertcat(satelite.state_vars(1:14,1)));
ee_fun = casadi.Function('ee',{opt.model.states(1:6+satelite.robot.n_q)},{R'*ee_fun(:,end)});

opt.costs.stage.parameters = opt.parameters.name(1:4);
opt.costs.stage.function   = @(x,u,varargin) sum((ee_fun(x(1:6+satelite.robot.n_q))-varargin{:}(4:6)).^2) + varargin{:}(7)*( x(18:22)'*x(18:22) + 1e2*harmf(x(1:22))'*harmf(x(1:22)) + (hsatf(x(1:22))-varargin{:}(8:10).*[0.5;0.5;0.5])'*(hsatf(x(1:22))-varargin{:}(8:10).*[0.5;0.5;0.5])    )*1e2 + 1e3*mu2(x) ;

opt.costs.general.parameters = opt.parameters.name(1:2);
opt.costs.general.function   = @(x,u,varargin) 1e3*(varargin{end}-varargin{end-1})'*(varargin{end}-varargin{end-1});

opt.constraints.states.upper  = vertcat( inf*ones(3,1),  inf*ones(3,1), 3600*2*pi/60*ones(3,1),  inf*ones(satelite.robot.n_q-3,1),  inf*ones(3,1),  0.09*ones(satelite.robot.n_q-3,1));
opt.constraints.states.lower  = vertcat(-inf*ones(3,1), -inf*ones(3,1), -3600*2*pi/60*ones(3,1), -inf*ones(satelite.robot.n_q-3,1), -inf*ones(3,1), -0.09*ones(satelite.robot.n_q-3,1));

opt.constraints.control.upper = vertcat(0.175*ones(3,1),50*ones(5,1));
opt.constraints.control.lower = -opt.constraints.control.upper;

opt.constraints.general.parameters  = opt.parameters.name(2);
opt.constraints.general.function{1} = @(x,varargin) (ee_fun(x(1:6+satelite.robot.n_q))-varargin{:}(1:3));
opt.constraints.general.elements{1} = 'end';
opt.constraints.general.type{1} = 'equality';

opt.constraints.parameters.name  = opt.parameters.name(2);
opt.constraints.parameters.upper =  vertcat(inf*ones(3,1));
opt.constraints.parameters.lower = -vertcat(inf*ones(3,1));

% Define inputs to optimization
opt.input.vector = opt.parameters.name([1 3 4]);
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
for k = 1:20000
    k

    link_positions(:,:,k) = full(R0s(:,:,k)'*satelite.kinematics.rL(R0s(:,:,k),vertcat(xsat(1:6+satelite.robot.n_q,k))));
    ref =  [2.806; -0.9745; -0.1926];

    if norm(link_positions(:,end,k)-ref)>=0.0005 && sw == 0
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

    mpc_input = vertcat(xsat(:,k),ref, auxi,sign(xsat(1:3,k)));
%     mpc_input = vertcat(xsat(:,k),ref, auxi,sign(pred_states(1:3,end)));
    sol = solver_mpc('x0', mpc_x0, 'lbx', args_mpc.lbx, 'ubx', args_mpc.ubx, ...
        'lbg', args_mpc.lbg, 'ubg', args_mpc.ubg, 'p', mpc_input);
    mpc_x0        = full(sol.x);
    torque(:,k) = full(sol.x(args_mpc.vars{3}));
    pred_states = reshape(full(sol.x(1:length(opt.model.states)*opt.N)),length(opt.model.states),opt.N);
    feval = full(opt.model.function(xsat(:,k),torque(:,k)));
    xsat(:,k+1) = xsat(:,k) + opt.dt*feval; %vertcat(xsat(satelite.idx.velocities,k),full(feval));
    
%     check_feas(solver_mpc.stats())


    tt(:,k) = full(dq0(R0s(:,:,k),xsat(1:22,k)));
    % update satelite quaternion, R0 according to next omega
    [q0s(:,k+1), R0s(:,:,k+1)]  = quaternion.integrate(tt(1:3,k),q0s(:,k),opt.dt);
    Rtest = Rf(xsat(1:3,k));
    xsat_full = vertcat(xsat(1:14,k),tt(:,k),xsat(15:22,k));
    mom_arm(:,k) = full(moment_arm(R0s(:,:,k),xsat_full));
    mom_sate(:,k) = full(moment_satelite(R0s(:,:,k),xsat_full));
    mom_wheels(:,k) = full(moment_wheels(R0s(:,:,k),xsat_full));
    itmt(:,k) = full(sol.x(end-5:end-3));
    for jj = 1:opt.N
        pred_mom_arm(:,jj) = full(harmf(pred_states(1:22,jj)));
    end
    tt2(:,k) = pred_mom_arm(:,3);

end

%%
close all

figure
% plot(itmt(1,:),'-.r')
hold on
plot(ref(1,:)*ones(1,length(squeeze(link_positions(1,end,:)))),'--r')
plot(squeeze(link_positions(1,end,:)),'-r')
% plot(itmt(2,:),'-.g')
plot(ref(2,:)*ones(1,length(squeeze(link_positions(1,end,:)))),'--g')
plot(squeeze(link_positions(2,end,:)),'-g')
% plot(itmt(3,:),'-.b')
plot(ref(3,:)*ones(1,length(squeeze(link_positions(1,end,:)))),'--b')
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
plot(torque(1,:),'r')
hold on
plot(torque(2,:),'g')
plot(torque(3,:),'b')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on
title('Reaction wheels torques')

figure
stairs(torque(4,:),'r')
hold on
stairs(torque(5,:),'g')
stairs(torque(6,:),'b')
stairs(torque(7,:),'c')
stairs(torque(8,:),'m')
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on
title('Manipulator joint torques')

% figure
% plot(xsat(23,:),'r')
% hold on
% plot(xsat(24,:),'g')
% plot(xsat(25,:),'b')
% title('Integral of h_{sat}')


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
plot(xsat(10,:),'r')
hold on
plot(xsat(11,:),'g')
plot(xsat(12,:),'b')
plot(xsat(13,:),'k')
plot(xsat(14,:),'m')
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