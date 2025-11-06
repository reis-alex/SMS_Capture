clear, clc, close all
import casadi.*

freefloating_path = [pwd '\floating_body.urdf'];
[freefloating,~] = urdf2robot_flex_visu(freefloating_path);
[freefloating_sim]= robot_slx_format(freefloating);

freefloating = SPART_casadi(freefloating_path);

% auxiliary skew symmetric matrix
SkewSym = @(x)[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

% for quaternion integration and conversions
quaternion = quaternion();

ddX_ff = freefloating.ddX();
kin_ff = freefloating.kinematics();

xff = [zeros(6,1); 
       1;0;1;
       1;0;0];                        

dt = 0.1;

% initialize satelite's state vector
theta_ff       = xff(1:3,1);
pos_ff          = xff(4:6,1);
omega0_ff      = xff(7:9,1);
r0dot_ff       = xff(10:12,1);

R0s = eye(3);
q0s = quaternion.rotm_to_quat(R0s); %initial quaternion, obtained directly from R0
theta_ff(:,1) = quaternion.quat_to_angles(q0s);
pos_f_CoM(:,1) = full(kin_ff.rLf(R0s,xff(4:6,1),0)); %just to check

for k = 1:100
    pos_ff          = xff(4:6,k);
    omega0_ff      = xff(7:9,k);
    r0dot_ff       = xff(10:12,k);

    feval = ddX_ff(R0s,pos_ff,omega0_ff,0,r0dot_ff,0,0);
    xff(:,k+1) = xff(:,k) + dt*vertcat(xff(7:end,k),full(feval));


    % update satelite quaternion, R0 according to next omega
    [q0s, R0s]  = quaternion.integrate(xff(7:9,k+1),q0s,dt);
    theta_ff(:,k+1)      = quaternion.quat_to_angles(q0s);
    pos_f_CoM(:,k+1) = full(kin_ff.rLf(R0s,xff(4:6,k+1),0)); %just to check
end


%% test jacobian
J        = freefloating.Jacob(1);
full([J.J0f(R0s,pos_f_CoM(:,end),0) J.Jmf(R0s,pos_ff(:,end),0)])