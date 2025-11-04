clear all;
close all;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Choice of robtos %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Robot from urdf Description %%%

[robot_post,robot_post_key] = urdf2robot_flex_visu('robot_post.urdf');
[robot_pre,robot_pre_key] = urdf2robot_flex_visu('robot_pre.urdf');
[target,target_key] = urdf2robot_flex_visu('target.urdf');


[robot_post_sim]= robot_slx_format(robot_post);
[robot_pre_sim]= robot_slx_format(robot_pre);
[target_sim]= robot_slx_format(target);

cmpt =1;
for i= 1:robot_post.n_links_joints
    if contains(robot_post.joints(i).name, 'Joint_6')
        id_EE(cmpt) = i;
        cmpt = cmpt+1;
    end
end
robot_post_sim.id_EE = id_EE;

cmpt =1;
for i= 1:robot_pre.n_links_joints
    if contains(robot_pre.joints(i).name, 'Joint_6')
        id_EE(cmpt) = i;
        cmpt = cmpt+1;
    end
end
robot_pre_sim.id_EE = id_EE;

cmpt =1;
for i= 1:target.n_links_joints
    if contains(target.joints(i).name, 'Joint_6')
        id_EE(cmpt) = i;
        cmpt = cmpt+1;
    end
end
target_sim.id_EE = id_EE;

target_sim.mtot =mass_tot(target);

cmpt =1;
for i= 1:robot_post.n_links_joints
    if contains(robot_post.links(i).name, 'Target')
        id_target(cmpt) = i;
        cmpt = cmpt+1;
    end
end
robot_post_sim.id_target = id_target(1);
robot_pre_sim.id_target = id_target(1);

cmpt =1;
for i= 1:target.n_links_joints
    if contains(target.links(i).name, 'Target')
        id_target(cmpt) = i;
        cmpt = cmpt+1;
    end
end
target_sim.id_target = id_target(1);


%--- System_dimensions ---%
dim_0 = 6;                                  % Base's number of DoFs 
dim_q = robot_post_sim.n_q;                      % Number of base+manipulator(s)'s actuators
dim_f = robot_post_sim.flexs.nb_mode;            % Number of flexible modes
robot_post_sim.dimensions.dim_r = 3;             % Number of reaction-wheels  
robot_post_sim.dimensions.dim_m = dim_q-robot_post_sim.dimensions.dim_r;   % Number of manipulator's actuators
robot_post_sim.dimensions.nb_manip = 1;          % Number of manipulators
dim_r = robot_post_sim.dimensions.dim_r;
robot_post_sim.dimensions.dim_m1 = 6;
dim_m1 = robot_post_sim.dimensions.dim_m1 ;
robot_pre_sim.dimensions = robot_post_sim.dimensions;
qr_init  = zeros(dim_r,1);

qv1_init = [0;0.524;-1.31;0;-0.785;0];
qv1_end = qv1_init+[pi/2;-pi/2;pi;0;pi;0]*0.01;

q_init = [qr_init;qv1_init+[pi/2;-pi/2;pi;0;pi;0]*0.2];
q_init = [qr_init;qv1_init];
q_end  = [qr_init;qv1_end];

qrdot_init  = 0*150*ones(dim_r,1);
qv1dot_init = zeros(dim_m1,1);
qdot_init = [qrdot_init; qv1dot_init];

qv1_init = [0;0.524;-1.31;0;-0.785;0]+[pi/2;-pi/2;pi;0;pi;0]*0.1;
qv1_end = qv1_init+[pi/2;-pi/2;pi;0;pi;0]*0.1;
q_init = [qr_init;qv1_init];
q_end  = [qr_init;qv1_end];

% d1 = display_robot(target,q_init);
d3 = display_robot(robot_post,q_init);
d5 = display_robot(robot_post,q_end);
% d4 = display_robot(robot_pre,q_init);


%%
%--- Initial conditions ---%

R0      = eye(3);                           % Base's intial rotation from Rsat to Rine
r0      = zeros(3,1);                       % Base's intial position in Rsat
q00     = [0;0;0;1;zeros(3,1)];             % Base's intial state=>[quaternion;r0]  
qm0     = q_init;                           % Actuators's initial value  
eta0    = zeros(dim_f,1);                   % Flexible modes' initial vibration
u00     = zeros(6,1);                       % Base's initial velocity in Rine
um0     = qdot_init;                         % Actuators' initial velocity in Rine
etadot0 = zeros(dim_f,1);                   % Flexible modes' initial velocity
qm = qm0; 
q0_dot = u00;
um = um0;
etadot = etadot0;
eta = eta0;
etaddot=zeros(dim_f,1);

% u0_t_0 = [0.001*ones(3,1); 0.001*ones(3,1)]*0.1;
u0_t_0 = rand(6,1)*1e-4;

E_inertia_max   = 0.05;
E_mass_max      = 0.05;
E_mode_max      = 0.05;
E_mass_target_max = 0.15;
E_inertia_target_max = 0.15;

robot_unc_pre        = get_robot_unc(robot_pre_sim,robot_pre,E_inertia_max,E_mass_max,E_mode_max,E_mass_target_max,E_inertia_target_max);
robot_unc_post       = get_robot_unc(robot_post_sim,robot_post,E_inertia_max,E_mass_max,E_mode_max,E_mass_target_max,E_inertia_target_max);



%%-- Control  --%%
wn    = 1;
ksi   = 0.7;
Kp   = wn*ksi*ones(dim_q,1);
Kv    = 2*ksi*wn*ones(dim_q,1);

robot_unc_pre.control.Kp = diag([Kp])/10;

robot_unc_pre.control.Kv = robot_unc_pre.control.Kp/100;

%%-- path planning %%
T_end           = 2000;
% T_end = 1500;
T_rest          = 500;
T_rest_init     = 0;
n_steps         = (T_rest_init+T_end+T_rest)/0.1;

[qm_ref, qmdot_ref,qmddot_ref, time] = path_planning(qv1_init.',qv1_end.', T_end, T_rest,T_rest_init, n_steps,'arm 1');

qm_ref_simu.time      = time;
qm_ref_simu.signals.dimensions     = dim_m1;
qm_ref_simu.signals.values         = [qm_ref.'].';

qmdot_ref_simu.time      = time;
qmdot_ref_simu.signals.dimensions     = dim_m1;
qmdot_ref_simu.signals.values         = [qmdot_ref.'].';



% robot = importrobot('robot_pre.urdf');
% robot.DataFormat = 'row';
% for i = 1:length(qm_ref)
%     show(robot,horzcat(qm_ref(i,:),zeros(1,3)));
%     drawnow
%     axis auto
% end
