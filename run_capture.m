addpath(genpath([pwd '\Toolbox_SMS_sym']));
clear all, clc
% close all
import casadi.*

% load robot through URDFs: Satellite + manipulator before capture
robot_path_pre = [pwd '\sms_simple_iftomm.urdf'];
[robot_pre,~] = urdf2robot_flex_visu(robot_path_pre);
[robot_pre_sim]= robot_slx_format(robot_pre);

% auxiliary skew symmetric matrix
SkewSym = @(x)[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

% for quaternion integration and conversions
quaternion = quaternion();

%% System dimensions
xsat = vertcat(zeros(6,1), ...                  % null initial displacements
    [pi/2, pi/2, pi/2, pi/8, pi/2, pi/4]', ...  % first 3: q_wheels, last 3: q_manip
    zeros(12,1));                               % null velocities

% initialize satelite's state vector
theta_sat       = xsat(1:3,1);
pos_sat         = xsat(4:6,1);
q_sat           = xsat(7:12,1);
omega0_sat      = xsat(13:15,1);
r0dot_sat       = xsat(16:18,1);
qrdots          = xsat(19:21,1);
qmdots          = xsat(22:24,1);
R0s = eye(3);
q0s = quaternion.rotm_to_quat(R0s); %initial quaternion, obtained directly from R0

% generate target dynamics
target_params.m = 1;
target_params.I = diag([1,2,3]);
target = free_body_euler(target_params);

% generate robot dynamics (using SPART)
robot_pre = SPART_casadi(robot_path_pre);

% Acceleration functions
ddX_target = target.ddX_ine();
ddX_robot = robot_pre.ddX();

Ht = target.ine_H();
Hs = robot_pre.H(); Hs = Hs.Hf;

kin_sat = robot_pre.kinematics();
link_positions = kin_sat.rLf(R0s,pos_sat,q_sat);
pos_t = full(link_positions(:,end));

% initialize target's state vector
R0t = eye(3);
q0t = quaternion.rotm_to_quat(R0t);
theta_t(:,1) = quaternion.quat_to_angles(q0t);
rdot_t(:,1)  = [0;0;0];
omega_t(:,1) = [0;0;1];
x_target = vertcat(pos_t, omega_t(:,1), rdot_t(:,1));

% simulation parameters

N = 500;
is_captured(1:N) = 0;
is_captured(100:end) = 1;
dt = 0.1;

% simulation loop

for k = 1:500
    if ~is_captured(:,k)
        % simulate free target
        xddot_target = full(ddX_target(R0t,omega_t(:,k),rdot_t(:,k),zeros(6,1)));
        x_target(:,k+1) = x_target(:,k) + dt*vertcat(x_target(7:9,k),xddot_target);

        % update target quaternion, R0 according to next omega
        [q0t, R0t]  = quaternion.integrate(x_target(4:6,k+1),q0t,dt);
        theta_t(:,k+1) = quaternion.quat_to_angles(q0t);

        pos_t(:,k+1)   = x_target(1:3,k+1);
        omega_t(:,k+1) = x_target(4:6,k+1);
        rdot_t(:,k+1)  = x_target(7:9,k+1);
    
        % simulate free robot
        pos_sat        = xsat(4:6,k);
        q_sat          = xsat(7:12,k);
        omega0_sat     = xsat(13:15,k);
        r0dot_sat      = xsat(16:18,k);
        qdot_sat       = xsat(19:24,k);

        torque_r(:,k) = zeros(3,1);
        torque_m(:,k) = zeros(3,1);
        feval = ddX_robot(R0s,pos_sat,omega0_sat,q_sat,r0dot_sat,qdot_sat,vertcat(torque_r(:,k),torque_m(:,k)));
        xsat(:,k+1) = xsat(:,k) + dt*vertcat(xsat(13:end,k),full(feval));

        % update satelite quaternion, R0 according to next omega
        [q0s, R0s]  = quaternion.integrate(xsat(13:15,k+1),q0s,dt);
        theta_sat(:,k)      = quaternion.quat_to_angles(q0s);

    else
        % Recall: dX_captured is vertcat(q0dot,qdot,qtdot)
        [dX_captured, ~] = freefloating_is_captured(xsat(:,k), R0s, x_target(:,k), R0t, robot_pre, robot_pre_sim, target, zeros(6,1));
        
        % update satelite's velocities, quaternion, and then integrate
        xsat(13:end,k+1) = full(dX_captured(1:12));
        
        %update quaternion and R0
        [q0s, R0s]  = quaternion.integrate(xsat(13:15,k+1),q0s,dt);

        % update satelite's positions
        xsat(1:12,k+1) = xsat(1:12,k) + dt*xsat(13:end,k);
        theta_sat(:,k)      = quaternion.quat_to_angles(q0s);

        % update target velocities, quaternion
        x_target(4:9,k+1) = full(dX_captured(13:end));
        [q0t, R0t]  = quaternion.integrate(x_target(4:6,k+1),q0t,dt);

        % integrate 
        x_target(1:3,k+1) = x_target(1:3,k) + dt*x_target(7:9,k);
        theta_t(:,k+1) = quaternion.quat_to_angles(q0t);

        % this is only for debugging, trying to see if shit works
        link_positions = kin_sat.rLf(R0s,pos_sat,q_sat);
        pos_ee = full(link_positions(:,end));
        [pos_ee x_target(1:3,end)]
                
        Hs_v = full(Hs(R0s,pos_sat,q_sat));
        Ht_v = full(Ht(R0t));
        [Hs_v(1:6,1:6)*xsat(13:18,k+1) Hs_v(7:end,7:end)*xsat(19:end,k+1) Ht_v*x_target(4:end,k+1)]
    end
end