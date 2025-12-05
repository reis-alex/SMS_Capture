
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
plot(torque(1,:));
hold on;
plot(torque(2,:))
plot(torque(3,:))
xlabel('Time [s]')
ylabel('Torque [Nm]')
grid on
title('Reaction wheels torques')

figure
stairs(torque(1,:),'r')
hold on
stairs(torque(2,:),'g')
stairs(torque(3,:),'b')
stairs(torque(4,:),'m')
stairs(torque(5,:),'c')
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
