addpath(genpath([pwd '\Toolbox_SMS_sym']));
clear all, clc
close all
import casadi.*

% load robot through URDFs, generate model
robot_path_pre = [pwd '\sms_simple_iftomm.urdf'];
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
    zeros(satelite.robot.n_q-3,1), pi/2,-0.1569,0, ...  % first 3: q_wheels, last 3: q_manip
    zeros(length(satelite.idx.velocities),1));                               % null velocities


% initialize satelite's state vector
R0s = eye(3);
q0s = quaternion.rotm_to_quat(R0s); %initial quaternion, obtained directly from R0


Ht = target.dynamics.H;
Hs = satelite.dynamics.H;
Cs = satelite.dynamics.C;

link_positions(:,:,1) = full(satelite.kinematics.rL(R0s,xsat));
pos_t = full(link_positions(:,end,1))
%%

% initialize target's state vector
R0t(:,:,1) = eye(3);
q0t(:,1) = quaternion.rotm_to_quat(R0t(:,:,1));
theta_t(:,1) = quaternion.quat_to_angles(q0t(:,1));
rdot_t(:,1)  = [0.1;0;0];
omega_t(:,1) = [0;0.0;0];
x_target = vertcat(theta_t, pos_t, omega_t(:,1), rdot_t(:,1));
theta_ff(:,1) = quaternion.quat_to_angles(q0t(:,1));
pos_f_CoM(:,1) = full(target.kinematics.rL(R0t(:,:,1),x_target)); %just to check

% simulation parameters

N = 500;
dt = 0.01;
total_moment = [];
kcapt = 1;


%% Create NDI: reaction wheels
Kp1 = 1.5*eye(3);
Kv1 = Kp1*0.9;
K1  = horzcat(Kv1,Kp1);

Hs2 = Hs(R0s,satelite.state_vars);
Hwr = Hs2(1:3,7:9);
Hr  = Hs2(7:9,7:9);
Hw  = Hs2(1:3,1:3);

Cs2 = Cs(R0s,satelite.state_vars);
Cr0 = Cs2(1:3,4:6); % 3x6
Cr  = Cs2(1:3,1:3);
C0r = Cr0'; % 6x3
C0 = Cs2(4:6,4:6);

iH0r  = [inv(Hwr)];
Hstar = Hwr' - Hr*iH0r*Hw;
Cstar = Cr0 - Hr*iH0r*C0 + Hr*iH0r*C0r*iH0r*H0(1:3,1:3) - Cr*iH0r*H0(1:3,1:3);

theta0   = satelite.state_vars(1:3);
q0des    = SX.sym('q0des',3,1);
q0dotdes = SX.sym('q0dotdes',3,1);
wr = Hstar*K1*[q0dotdes-omega0;q0des-theta0(1:3)]+ Cstar(:,1:3)*theta0(1:3);

wr = Function('NDI_wheels',{R0,satelite.state_vars,q0des,q0dotdes},{wr});


%% simulation loop
            EE_idx = find(strcmp({satelite.robot.links.name},'Link_EE'));
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
        torque_r(:,k) = full(wr(R0s(:,:,k),xsat(:,k), zeros(3,1), zeros(3,1)));
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

            hcapt = compute_moment(satelite, target,R0s(:,:,k),R0t(:,:,k),xsat(:,k),x_target(:,k));

            Je  = full(Jee.Jmf(R0s(:,:,kcapt),xsat(:,kcapt)));
            J0  = full(Jee.J0f(R0s(:,:,kcapt),xsat(:,kcapt)));
            Jt = blkdiag(R0t(:,:,kcapt),eye(3));
            if satelite.robot.n_q>3
                xsat(satelite.idx.velocities,kcapt)  = full([P0*Hs_v(1:6,1:6) P0*Hs_v(1:6,7:end); J0 Je]\[zeros(6,1); +Jt*x_target(target.idx.velocities,kcapt)]);
            else
                xsat(satelite.idx.velocities,kcapt)  = full(pinv([P0*Hs_v(1:6,1:6) P0*Hs_v(1:6,7:end); J0 Je])*[zeros(6,1); +Jt*x_target(target.idx.velocities,kcapt)]);
            end
        end

        % Recall: dX_captured is vertcat(q0dot,qdot,qtdot)
        torque_r(:,k) = zeros(satelite.robot.n_q-3,1);
        torque_m(:,k) = zeros(3,1);
        [xsat_plus, xtarget_plus] = freefloating_is_captured(dt, xsat(:,k), R0s(:,:,k), x_target(:,k), R0t(:,:,k), satelite, robot_pre_sim, target, vertcat(torque_r(:,k),torque_m(:,k)), hcapt);
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