close all;clc;

% step = 300;
% display_robot_motion_simu(robot_post,out.logsout,step);
% 
log = out.logsout;
R0 = log.getElement('R0');
r0 = log.getElement('r0');
qm = log.getElement('q_simu');
rL_target = log.getElement('rL_target');
rL_EE = log.getElement('rL_EE');
rL_EE_ref = log.getElement('rL_EE_ref');
u_EE = log.getElement('u_EE');
u_target = log.getElement('u_t');
w_ref = log.getElement('w_ref');
w = log.getElement('w');
q_ref = log.getElement('q_ref');
qdot_ref = log.getElement('qdot_ref');
q_simu = log.getElement('q_simu');
qdot_simu = log.getElement('qdot_simu');

time = R0.Values.Time;

% colorn = "gem12";
% colororder(colorn)
colorn={'red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'black','red', 'green', 'blue', 'cyan'};


figure,
subplot(2,1,1);
for i=1:dim_m1
    plot(time,qdot_ref.Values.Data(:,i),'LineWidth',2 ,'Color',colorn{i},'DisplayName', '$\dot{q}_{ref}$');
    hold on;
    plot(time,qdot_simu.Values.Data(:,dim_r+i),'-.','LineWidth',2,'Color',colorn{i},'DisplayName', '$\dot{q}_{simu}$');
    hold on;
end
%plot(time,qdot_simu.Values.Data(:,dim_r+1:end),'-.','DisplayName', '$\dot{q}_{simu}$');
title('Comparison of manipualtor velocities path planning vs simu');
hl = legend('show');
set(hl, 'interpreter','tex'); grid on;
subplot(2,1,2);
for i=1:dim_m1
    plot(time,q_ref.Values.Data(:,i),'LineWidth',2 ,'Color',colorn{i},'DisplayName', '$q_{ref}$');
    hold on;
    plot(time,q_simu.Values.Data(:,dim_r+i),'-.','LineWidth',2 ,'Color',colorn{i},'DisplayName', '$q_{simu}$');
    hold on;
end
title('Comparison of manipualtor poses path planning vs simu');
hl = legend('show');
set(hl, 'interpreter','tex'); grid on;


figure,
subplot(2,1,1);
for i=1:3
    plot(time,w_ref.Values.Data(i,:),'LineWidth',2 ,'Color',colorn{i},'DisplayName', 'w_{ref}');
    hold on;
    plot(time,w.Values.Data(i,:),'-.','LineWidth',2,'Color',colorn{i},'DisplayName', '$w_1$');
    hold on;
end
title('w tracking task 1');
hl = legend('show');
set(hl, 'interpreter','tex'); grid on;
subplot(2,1,2);
for i=4:robot_unc_pre.n_q
    plot(time,w_ref.Values.Data(i,:),'LineWidth',2 ,'Color',colorn{i},'DisplayName', '$w_ref$');
    hold on;
    plot(time,w.Values.Data(i,:),'-.','LineWidth',2,'Color',colorn{i},'DisplayName', '$w_1$');
    hold on;
end
title('w tracking task 2');
hl = legend('show');
set(hl, 'interpreter','tex'); grid on;


figure,
plot(time,squeeze(rL_EE_ref.Values.Data),'LineWidth',2,'DisplayName', '$rL_{EE_{ref}}$')
hold on;
plot(time,rL_EE.Values.Data,'LineWidth',2,'Marker','.','DisplayName', '$rL_{EE}$')
title('End-Effector pose tracking');
hl = legend('show');
set(hl, 'interpreter','tex'); grid on;















%% UTILS

function data_out = get_data_in(data,name)
    for i=1:data.numElements
        data_tmp = data{i}.Values;
        if size(data_tmp.Name) == size(name)
            if data_tmp.Name == name
                data_out = data{i}.Values;
            end
        end
    end
end
