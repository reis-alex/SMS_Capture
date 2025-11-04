clc;close all;
nb_run =12;


% for i =8:nb_run
%     % u0_t_0 = rand(6,1)*1e-4;
%     u0_t_0 = [rand(3,1)*1e-4; zeros(3,1)];
% 
%     sim("Robust_task_NDI_capture.slx");
%     save(strcat("NDI_NDO_run_",int2str(i)),"ans");
%     clear ans;
% 
% end

%%


close all;
% step = 300;
% % % display_robot_motion_simu(robot_post,out.logsout,step);
% display_robot_motion_simu(robot_post,ans.logsout,step,t_start);

Fs = 25;

for i=1:nb_run
    load(strcat("NDI_NDO_run_",int2str(i)));
    run(i) = ans.logsout; 

    R0(i) = run(i).getElement('R0');
    r0(i) = run(i).getElement('r0');
    qm(i) = run(i).getElement('q_simu');
    rL_target(i) = run(i).getElement('rL_target');
    rL_EE(i) = run(i).getElement('rL_EE');
    rL_EE_ref(i) = run(i).getElement('rL_EE_ref');
    u_EE(i) = run(i).getElement('u_EE');
    u_target(i) = run(i).getElement('u_t');
    w_ref(i) = run(i).getElement('w_ref');
    w(i) = run(i).getElement('w');
    q_ref(i) = run(i).getElement('q_ref');
    qdot_ref(i) = run(i).getElement('qdot_ref');
    q_simu(i) = run(i).getElement('q_simu');
    qdot_simu(i) = run(i).getElement('qdot_simu');
    q0dot_simu(i)= run(i).getElement('q0dot_simu');
    Fe(i)= run(i).getElement('Fe');
    rcom(i) = run(i).getElement('r_com');

end

time = R0(1).Values.Time;



t_start = round(length(time)/time(end)*2000)+1;

figure,
for i =1:nb_run
    % subplot(2,1,2),
    % plot(time,squeeze(w_ref(i).Values.Data)-squeeze(w(i).Values.Data),DisplayName=strcat('run',int2str(i)), LineWidth=3);
    % hold on;
    subplot(2,1,1),
    plot(time,squeeze(w(i).Values.Data),DisplayName=strcat('run',int2str(i)), LineWidth=3);
    hold on;
end 
grid on;
%legend;
title('Control tracking performances', 'Interpreter', 'latex','FontSize',Fs);
xlabel('time (s)', 'Interpreter', 'latex','FontSize',Fs);
ylabel('(m/s) or (rad/s)', 'Interpreter', 'latex','FontSize',Fs);
ax = gca;
ax.LineWidth=2;
ax.FontSize = Fs;
xlim([0 3500]);

for i =1:nb_run
    subplot(2,1,2),
    plot(time,squeeze(w_ref(i).Values.Data)-squeeze(w(i).Values.Data),DisplayName=strcat('run',int2str(i)), LineWidth=3);
    hold on;

end 
grid on;
%legend;
title('Control error signal $\dot{\mathbf{\epsilon}}_c$', 'Interpreter', 'latex','FontSize',Fs);
xlabel('time (s)', 'Interpreter', 'latex','FontSize',Fs);
ylabel('(m/s) or (rad/s)', 'Interpreter', 'latex','FontSize',Fs);
ax = gca;
ax.LineWidth=2;
ax.FontSize = Fs;
xlim([0 3500]);


figure,
for i =1:nb_run
    plot(time, sqrt(sum(u_target(i).Values.Data.^2,2)),DisplayName=strcat('run',int2str(i)), LineWidth=3);
    hold on;

end 
grid on;
legend;
title('Evolution of the norm of the target velocities', 'interpreter','tex','FontSize',Fs);
ax = gca;
ax.LineWidth=2;
ax.FontSize = Fs;
xlabel('time (s)', 'Interpreter', 'latex','FontSize',Fs);
ylabel('$\| \mathbf{t}_{t} \|$', 'Interpreter', 'latex','FontSize',Fs);
xlim([0 3500]);

% figure,
% for i =1:nb_run
%     plot(time, sqrt(sum(q0dot_simu(i).Values.Data.^2,2)),DisplayName=strcat('run',int2str(i)), LineWidth=3);
%     hold on;
% 
% end 
% grid on;
% legend;
% title('norm u_{0}', 'interpreter','tex');

figure,
for i =1:nb_run
    plot(time, sqrt(sum(q0dot_simu(i).Values.Data(:,1:3).^2,2)),DisplayName=strcat('run',int2str(i)), LineWidth=3);
    hold on;

end 
grid on;
legend;
%title('norm \omega_{0}', 'interpreter','tex');
%title('Evolution of \omega_{0} norm', 'interpreter','tex');
title('Evolution of the norm of the base angular velocity', 'interpreter','tex');
ax = gca;
ax.LineWidth=2;
ax.FontSize = Fs;
xlabel('time (s)', 'Interpreter', 'latex','FontSize',Fs);
ylabel('$\| \omega_{0} \|$', 'Interpreter', 'latex','FontSize',Fs);
xlim([0 3500]);



figure,
for i=1:nb_run
    plot3(rcom(i).Values.Data(:,1),rcom(i).Values.Data(:,2),rcom(i).Values.Data(:,3),'LineWidth',3,'LineStyle','-','DisplayName','r_0');
    hold on;
    plot3(r0(i).Values.Data(:,1),r0(i).Values.Data(:,2),r0(i).Values.Data(:,3),'LineWidth',3,'LineStyle',':','DisplayName','r_0');
    hold on;
    plot3(rL_EE(i).Values.Data(:,1),rL_EE(i).Values.Data(:,2),rL_EE(i).Values.Data(:,3),'LineWidth',3,'LineStyle','-.','DisplayName','rL_EE');
end
grid on;
%legend;
title('Evolution of base, end-effector and CoM pose in Rine','FontSize',Fs);
xlabel('X (m)', 'Interpreter', 'latex','FontSize',Fs);
ylabel('Y (m)', 'Interpreter', 'latex','FontSize',Fs);
zlabel('Z (m)', 'Interpreter', 'latex','FontSize',Fs);


figure,
sgtitle('Relative displacement of base/End-effectors during the full motion (left) and the post-capture (right)', 'interpreter','tex','FontSize',Fs);
subplot(1,2,1),
for i=1:nb_run
    plot(time,sqrt(sum(r0(i).Values.Data.^2,2)),'LineWidth',3,'LineStyle','-','DisplayName',strcat('run',int2str(i)));
    hold on;
end
hold on;
for i=1:nb_run
    hold on;
    plot(time,sqrt(sum(rL_EE(i).Values.Data.^2,2))-norm(rL_EE(i).Values.Data(1,:)),'LineWidth',3,'LineStyle',':','DisplayName',strcat('run',int2str(i)));
end
l =legend; 
l.FontSize = Fs/2;
grid on;
% title('Evolution of the norm of the base angular velocity', 'interpreter','tex');

ax = gca;
ax.LineWidth=2;
ax.FontSize = Fs;
xlabel('time (s)', 'Interpreter', 'latex','FontSize',Fs);
ylabel('$\| r_{\mathcal{S}_{0}} \|$ full line (m), $\| r_{L_{EE}} \|$ dot line (m) ', 'Interpreter', 'latex','FontSize',Fs);
xlim([0 3500]);

subplot(1,2,2),
for i=1:nb_run
    plot(time(t_start:end),sqrt(sum(r0(i).Values.Data(t_start:end,:,:).^2,2)),'LineWidth',3,'LineStyle','-','DisplayName',strcat('run',int2str(i)));
    hold on;
end
hold on;
for i=1:nb_run
    hold on;
    plot(time(t_start:end),sqrt(sum(rL_EE(i).Values.Data(t_start:end,:,:).^2,2))-norm(rL_EE(i).Values.Data(t_start,:)),'LineWidth',3,'LineStyle',':','DisplayName',strcat('run',int2str(i)));
end
l =legend; 
l.FontSize = Fs/2;
grid on;
% title('Relative displacement of base/End-effectors', 'interpreter','tex','FontSize',Fs/2);
ax = gca;
ax.LineWidth=2;
ax.FontSize = Fs;
xlabel('time (s)', 'Interpreter', 'latex','FontSize',Fs);
ylabel('$\| r_{\mathcal{S}_{0}} \|$ full line (m), $\| r_{L_{EE}} \|$ dot line (m) ', 'Interpreter', 'latex','FontSize',Fs);
xlim([2000 3500]);


% 
% % load robot_post_capture.mat
% % load robot_pre_capture.mat
% % t_start = round(length(time)/time(end)*2500)+1;
% step = round((length(time)-t_start)/1);
% display_robot_motion_simu(robot_post,run(1),step,t_start);
% 
% display_robot_motion_simu(robot_pre,run(1),length(time)-1,1);

