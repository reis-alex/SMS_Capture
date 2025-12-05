addpath(genpath([pwd '\Toolbox_SMS_sym']));
clear all, clc
close all
import casadi.*

% load robot through URDFs, generate model
robot_path_pre = [pwd '\robot_isparo.urdf'];
[satelite,~] = urdf2robot_flex_visu(robot_path_pre);
[robot_pre_sim]= robot_slx_format(satelite);
satelite = SPART_casadi(robot_path_pre);
%%
ref =  [2.806; -0.9745; -0.1926];

load('accelerations.mat')
load('desired.mat')
load('EE.mat')
load('inter_ref.mat')
load('mom_arm.mat')
load('mom_sat.mat')
load('mom_wheel.mat')
load('NDI.mat')
load('states.mat')
load('tau_m.mat')

close all


figure
plot(itmt(1,:),'-.r')
hold on
plot(ref(1,:)*ones(1,length(squeeze(ee_pos(1,end,:)))),'--r')
plot(ee_pos(1,:),'-r')
plot(itmt(2,:),'-.g')
plot(ref(2,:)*ones(1,length(squeeze(ee_pos(1,end,:)))),'--g')
plot(ee_pos(2,:),'-g')
plot(itmt(3,:),'-.b')
plot(ref(3,:)*ones(1,length(squeeze(ee_pos(1,end,:)))),'--b')
plot(ee_pos(3,:),'-b')
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
stairs(torque_arm(4,:),'m')
stairs(torque_arm(5,:),'c')
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

for i=1:length(xsat(1,:))-1
   hdot(:,i) = satelite.robot.links(1).inertia*accelerations(1:3,i);
end
figure
plot(hdot(1,:),'r')
hold on
plot(hdot(2,:),'g')
plot(hdot(3,:),'b')