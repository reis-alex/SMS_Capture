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
clear x_target x_sat theta_t pos_t omega_t rdot_t
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
pos_t = full(link_positions(:,end,1))

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

r = - C0 + mtimes(G, Cq) - mtimes(S, virtual_control);

% NDI for wheels
G_r = G(:, 1:nr);
lambda = 1e-8;
tau_r = mtimes((G_r' * G_r + lambda*eye(nr)) \ (G_r' ), r);

NDI_wheels = Function('NDI_fun', {satelite.R0, satelite.state_vars, virtual_control}, {tau_r}, {'R0','x','v_des'}, {'tau_r'});

% NDI for manipulator
Kp2 = 1.5*eye(5);
Kv2 = Kp2*0.1;
K2  = horzcat(Kv2,Kp2);

virtual_control = casadi.SX.sym('v',8,1);
H00_reg_inv = inv(H00 + eps_reg*eye(6));
S = Hqq-mtimes(H0q',mtimes(H00_reg_inv,H0q));

tau_m = Cq-H0q'*H00_reg_inv*C0 + S*virtual_control;
tau_m = tau_m(4:end);
NDI_arm = Function('NDI_fun', {satelite.R0, satelite.state_vars, virtual_control}, {tau_m}, {'R0','x','v_des'}, {'tau_r'});

%% simulation loop
% simulation parameters
dt = 0.01;
total_moment = [];
feval = zeros(length(satelite.state_vars)/2,1);
for k = 1:20000
    k
    link_positions(:,:,k) = full(satelite.kinematics.rL(R0s(:,:,k),xsat(:,k)));

    % simulate free robot
    vr = -vertcat(K1*vertcat(xsat(1:3,k)-[0;0;0.1],xsat(satelite.idx.omega0,k)-[0;0;0]),zeros(3,1));
    torque_wheels(:,k) = vertcat(full(NDI_wheels(R0s(:,:,k),xsat(:,k),vr)), zeros(5,1));

    vr = -vertcat(zeros(3,1), K2*vertcat(xsat(satelite.idx.q(4:end),k)-[0;0;0;0;0.1],xsat(satelite.idx.qdot(4:end),k)-zeros(5,1)));
%     torque_arm(:,k) = zeros(8,1);
    torque_arm(:,k) = vertcat(zeros(3,1),full(NDI_arm(R0s(:,:,k),xsat(:,k),vr)));

    torque(:,k) = torque_wheels(:,k)+torque_arm(:,k);
    feval = satelite.dynamics.ddX(R0s(:,:,k),xsat(:,k),torque(:,k));
    xsat(:,k+1) = xsat(:,k) + dt*vertcat(xsat(satelite.idx.velocities,k),full(feval));

    % update satelite quaternion, R0 according to next omega
    [q0s(:,k+1), R0s(:,:,k+1)]  = quaternion.integrate(xsat(satelite.idx.omega0,k),q0s(:,k),dt);

end

%% 
colors = {'r','g','b','k','c','m'};

figure(1)
plot(xsat(1,:),'r')
hold on
plot(xsat(2,:),'g')
plot(xsat(3,:),'b')
title('Base orientation')


figure(2)
plot(xsat(4,:),'r')
hold on
plot(xsat(5,:),'g')
plot(xsat(6,:),'b')
title('Base position (r0)')

figure(3)
hold on
for i = satelite.idx.q
    plot(xsat(i,:),colors{i+1-satelite.idx.q(1)})
end
title('Servicer joint velocities (reaction wheels (--), arm joints (-))')
title('Servicer joints (reaction wheels (--), arm joints (-))')

figure(4)
hold on
for i = satelite.idx.velocities(1:3)
    plot(xsat(i,:),colors{i+1-satelite.idx.velocities(1)})
end
title('Base angular velocity')


figure(5)
hold on
for i = satelite.idx.velocities(4:6)
    plot(xsat(i,:),colors{i+1-satelite.idx.velocities(4)})
end
title('Base linear velocity (\dot(r)_0)')

figure(6)
hold on
for i = satelite.idx.qdot
    plot(xsat(i,:),colors{i+1-satelite.idx.qdot(1)})
end
title('Servicer joint velocities (reaction wheels (--), arm joints (-))')

figure(7)
plot(x_target(1,:),'r')
hold on
plot(x_target(2,:),'g')
plot(x_target(3,:),'b')
title('Target orientation')

figure(8)
plot(x_target(4,:),'r')
hold on
plot(x_target(5,:),'g')
plot(x_target(6,:),'b')
title('Target position')

figure(9)
plot(x_target(7,:),'r')
hold on
plot(x_target(8,:),'g')
plot(x_target(9,:),'b')
title('Target angular velocity')

figure(10)
plot(x_target(10,:),'r')
hold on
plot(x_target(11,:),'g')
plot(x_target(12,:),'b')
title('Target linear velocity')

EE_positions = squeeze(link_positions(:,end,:));
figure(11)
hold on
plot(EE_positions(1,:),'r')
plot(EE_positions(2,:),'g')
plot(EE_positions(3,:),'b')
plot(x_target(4,:),'--r')
plot(x_target(5,:),'--g')
plot(x_target(6,:),'--b')

%% Fun plot
figure
for i = 1:910
plot3(link_positions(1,:,i),link_positions(2,:,i),link_positions(3,:,i),'-o')
axis([-0.1 0.5 -1 1 -1 2.5])
drawnow
pause(0.01)
end

%% Auxiliary functions


function m = compute_moment(satelite, target, R0s,R0t,xsat,x_target)
    Ht_v = full(target.dynamics.H(R0t,x_target));
    Hs_v = full(satelite.dynamics.H(R0s,xsat));
    m = full(blkdiag(R0s,eye(3))*(Hs_v(1:6,1:6)*xsat(satelite.idx.omega0(1):satelite.idx.r0dot(end)) + Hs_v(1:6,7:end)*xsat(satelite.idx.qdot)) + blkdiag(R0t,eye(3))*Ht_v*x_target(target.idx.velocities));
end