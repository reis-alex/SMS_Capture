addpath(genpath([pwd '\Toolbox_SMS_sym']));
clear all, clc
close all
import casadi.*

% load robot through URDFs, generate model
robot_path_pre = [pwd '\robot_isparo.urdf'];
[satelite,~] = urdf2robot_flex_visu(robot_path_pre);
[robot_pre_sim]= robot_slx_format(satelite);
satelite = SPART_casadi(robot_path_pre);

% load target through URDF, generate model
freefloating_path = [pwd '\floating_body.urdf'];
[target,~] = urdf2robot_flex_visu(freefloating_path);
[freefloating_sim]= robot_slx_format(target);
target = SPART_casadi(freefloating_path);

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


Ht = target.dynamics.H;
Hs = satelite.dynamics.H;

link_positions(:,:,1) = full(satelite.kinematics.rL(R0s,xsat));
pos_t = full(link_positions(:,end,1))
%%

% initialize target's state vector
R0t(:,:,1) = eye(3);
q0t(:,1) = quaternion.rotm_to_quat(R0t(:,:,1));
theta_t(:,1) = quaternion.quat_to_angles(q0t(:,1));
rdot_t(:,1)  = [0;0;0];
omega_t(:,1) = [0;0.0;0.1];
x_target = vertcat(theta_t, pos_t, omega_t(:,1), rdot_t(:,1));
theta_ff(:,1) = quaternion.quat_to_angles(q0t(:,1));
pos_f_CoM(:,1) = full(target.kinematics.rL(R0t(:,:,1),x_target)); %just to check

% simulation parameters

N = 500;
dt = 0.1;
total_moment = [];
kcapt = 1;
% simulation loop
            EE_idx = find(strcmp({satelite.robot.links.name},'Link_6_'));
Jee = satelite.Jacob(EE_idx);

for k = 1:20000
    k
    link_positions(:,:,k) = full(satelite.kinematics.rL(R0s(:,:,k),xsat(:,k)));

    Je  = full(Jee.Jmf(R0s(:,:,k),xsat(:,k)));
    J0  = full(Jee.J0f(R0s(:,:,k),xsat(:,k)));
    Jt = blkdiag(R0t(:,:,k),eye(3));
    total_velocities(:,k) = vertcat(J0*xsat(satelite.idx.omega0(1):satelite.idx.r0dot(end),k)+Je*xsat(satelite.idx.qdot,k), ...
                                    Jt*x_target(target.idx.velocities,k));
    if k<kcapt
        % simulate free target
        feval = target.dynamics.ddX(R0t(:,:,k),x_target(:,k),0);
        x_target(:,k+1) = x_target(:,k) + dt*vertcat(x_target(7:end,k),full(feval));


        % update target quaternion, R0 according to next omega
        [q0t(:,k+1), R0t(:,:,k+1)]       = quaternion.integrate(x_target(target.idx.omega0,k+1),q0t(:,k),dt);
    
        % simulate free robot
        torque_r(:,k) = zeros(satelite.robot.n_q-3,1);
        torque_m(:,k) = zeros(satelite.robot.n_q,1);
        feval = satelite.dynamics.ddX(R0s(:,:,k),xsat(:,k),vertcat(torque_r(:,k),torque_m(:,k)));
        xsat(:,k+1) = xsat(:,k) + dt*vertcat(xsat(satelite.idx.velocities,k),full(feval));

        % update satelite quaternion, R0 according to next omega
        [q0s(:,k+1), R0s(:,:,k+1)]  = quaternion.integrate(xsat(satelite.idx.omega0,k),q0s(:,k),dt);
%         theta_sat(:,k)      = quaternion.quat_to_angles(q0s);

        total_moment = horzcat(total_moment, compute_moment(satelite, target,R0s(:,:,k),R0t(:,:,k),xsat(:,k),x_target(:,k)));


    else
        if k == kcapt
            Ht_v = full(target.dynamics.H(R0t(:,:,kcapt),x_target(:,kcapt)));
            Hs_v = satelite.dynamics.H(R0s(:,:,kcapt),xsat(:,kcapt));
            P0 = blkdiag(R0s(:,:,kcapt),eye(3));
            Pt = blkdiag(R0t(:,:,kcapt),eye(3));

            hcapt = compute_moment(satelite, target,R0s(:,:,kcapt),R0t(:,:,kcapt),xsat(:,kcapt),x_target(:,kcapt));

            Je  = full(Jee.Jmf(R0s(:,:,kcapt),xsat(:,kcapt)));
            J0  = full(Jee.J0f(R0s(:,:,kcapt),xsat(:,kcapt)));
            Jt = blkdiag(R0t(:,:,kcapt),eye(3));
%             if satelite.robot.n_q>3
%                 xsat(satelite.idx.velocities,kcapt)  = full([P0*Hs_v(1:6,1:6) P0*Hs_v(1:6,7:end); J0 Je]\[zeros(6,1); +Jt*x_target(target.idx.velocities,kcapt)]);
%             else
%                 xsat(satelite.idx.velocities,kcapt)  = full(pinv([P0*Hs_v(1:6,1:6) P0*Hs_v(1:6,7:end); J0 Je])*[zeros(6,1); +Jt*x_target(target.idx.velocities,kcapt)]);
%             end
        end

        % Recall: dX_captured is vertcat(q0dot,qdot,qtdot)
        torque_r(:,k) = zeros(satelite.robot.n_q-3,1);
        torque_m(:,k) = zeros(3,1);
        [xsat_plus, xtarget_plus] = freefloating_is_captured(dt, xsat(:,k), R0s(:,:,k), x_target(:,k), R0t(:,:,k), satelite, target, vertcat(torque_r(:,k),torque_m(:,k)), hcapt);
        xsat(:,k+1) = xsat_plus;
        x_target(:,k+1) = xtarget_plus;
        
        %update quaternion and R0
        [q0s(:,k+1), R0s(:,:,k+1)]  = quaternion.integrate(xsat(satelite.idx.omega0,k+1),q0s(:,k),dt);
        [q0t(:,k+1), R0t(:,:,k+1)]  = quaternion.integrate(x_target(target.idx.omega0,k+1),q0t(:,k),dt);

        % compute total moment at k
        total_moment = horzcat(total_moment, compute_moment(satelite, target, R0s(:,:,k),R0t(:,:,k),xsat(:,k),x_target(:,k)));

    end
end

%% 
colors = {'r','g','b','k','c','m'};

figure()
plot(xsat(1,:),'r')
hold on
plot(xsat(2,:),'g')
plot(xsat(3,:),'b')
title('Base orientation')

figure()
plot(xsat(4,:),'r')
hold on
plot(xsat(5,:),'g')
plot(xsat(6,:),'b')
title('Base position (r0)')

figure()
hold on
for i = 1:3
    plot(xsat(satelite.idx.r0(end)+i,:),colors{i})
end
title('Servicer joint  (reaction wheels)')

figure()
hold on
for i = satelite.idx.q(4:end)
    plot(xsat(i,:),colors{i+1-satelite.idx.q(4)})
end
title('Servicer joints (arm joints)')

figure()
hold on
for i = satelite.idx.velocities(1:3)
    plot(xsat(i,:),colors{i+1-satelite.idx.velocities(1)})
end
title('Base angular velocity')


figure()
hold on
for i = satelite.idx.velocities(4:6)
    plot(xsat(i,:),colors{i+1-satelite.idx.velocities(4)})
end
title('Base linear velocity (\dot(r)_0)')

figure()
hold on
for i = 1:3
    plot(xsat(satelite.idx.r0dot(end)+i,:),colors{i})
end
title('Servicer joint velocities (reaction wheels)')

figure()
hold on
for i = satelite.idx.qdot(4:end)
    plot(xsat(i,:),colors{i+1-satelite.idx.qdot(4)})
end
title('Servicer joints velocities (arm joints)')

figure()
plot(x_target(1,:),'r')
hold on
plot(x_target(2,:),'g')
plot(x_target(3,:),'b')
title('Target orientation')

figure()
plot(x_target(4,:),'r')
hold on
plot(x_target(5,:),'g')
plot(x_target(6,:),'b')
title('Target position')

figure()
plot(x_target(7,:),'r')
hold on
plot(x_target(8,:),'g')
plot(x_target(9,:),'b')
title('Target angular velocity')

figure()
plot(x_target(10,:),'r')
hold on
plot(x_target(11,:),'g')
plot(x_target(12,:),'b')
title('Target linear velocity')

EE_positions = squeeze(link_positions(:,end,:));
figure()
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
    grid on
    axis([-5 5 -5 5 -5 5])
    drawnow
    pause(0.01)
end

%% Auxiliary functions


function m = compute_moment(satelite, target, R0s,R0t,xsat,x_target)
    Ht_v = full(target.dynamics.H(R0t,x_target));
    Hs_v = full(satelite.dynamics.H(R0s,xsat));
    m = full(blkdiag(R0s,eye(3))*(Hs_v(1:6,1:6)*xsat(satelite.idx.omega0(1):satelite.idx.r0dot(end)) + Hs_v(1:6,7:end)*xsat(satelite.idx.qdot)) + blkdiag(R0t,eye(3))*Ht_v*x_target(target.idx.velocities));
end