function [dX, ddX] = freefloating_is_captured(xsat, R0s, xtarget, R0t, robot_pre, robot_pre_sim, target, tau)

% load state vector for satelite and target
r0_sat         = xsat(4:6);
q_sat          = xsat(7:12);
omega0_sat     = xsat(13:15);
r0dot_sat      = xsat(16:18);
qdot_sat       = xsat(19:24);

omega_target = xtarget(4:6);
rdot_target  =  xtarget(7:9,1);

% get matrices from the satelites dynamics (from SPART)
H = robot_pre.H();
H = full(H.Hf(R0s,r0_sat,q_sat));
Hq = H(7:12,7:12);
H0q = H(1:6,7:12);
H0 = H(1:6,1:6);

C = robot_pre.C();
C = full(C.Cf(R0s,r0_sat,q_sat,omega0_sat,r0dot_sat,qdot_sat))*[omega0_sat;r0dot_sat;qdot_sat];
C0 = C(1:6);
Cq = C(7:end);

% get matrices from the target dynamics (from free_body_euler)
Ht_ine  = target.ine_H();
Ht      = Ht_ine(R0t);
Ht_dot  = target.ine_H_dot();
H_bar   = blkdiag(H,Ht);

Ct_t    = target.ine_C();
Ct      = full(Ct_t(R0t,omega_target)*[omega_target; rdot_target]);

% target jacobian
Jt = target.J_point_inertial(zeros(3,1)); %zeros bc r_b is the distance to the CoM

% Jt_inv = inv(Jt(R0t)); 
Jt_inv = blkdiag(R0t',eye(3)); % after discussion with Mathieu (Jt computed wrt target's CoM, numerically better)

% end effector jacobian (from SPART)
EE_idx                  = find(strcmp({robot_pre_sim.links.name},'Link_EE'));
robot_kinematics        = robot_pre.kinematics(EE_idx);
robot_diffkinematics    = robot_pre.diffkinematics();

Jee = robot_pre.Jacob(EE_idx);
Je  = Jee.Jmf(R0s,r0_sat,q_sat);
J0  = Jee.J0f(R0s,r0_sat,q_sat);

Jee_dot = robot_pre.Jacobdot(EE_idx);

%% Build constrained, reduced dynamics (homogeneous with Overleaf)
Bcim = [eye(6) zeros(6);
        zeros(6) eye(6);
        -Jt_inv*J0 -Jt_inv*Je];


% get twist of the end effector rL(:,end) already in the inertial frame
% (from SPART)
robot_twist = robot_diffkinematics.tLf(R0s,r0_sat,omega0_sat,q_sat,r0dot_sat,qdot_sat);

% compute Je_dot and J0_dot from robot_pre
Je_dot = Jee_dot.Jmdotf(R0s,omega0_sat,r0_sat,q_sat,r0dot_sat,qdot_sat,robot_twist(:,end));
J0_dot = Jee_dot.J0dotf(R0s,omega0_sat,r0_sat,q_sat,r0dot_sat,robot_diffkinematics.t0f(R0s,omega0_sat,r0dot_sat));

% calc d(Jt_inv)/dt
Jt_dot = target.Jacobdot(zeros(3,1)); %zeros bc r_b is the distance to the CoM
Jt_inv_dot = -Jt_inv*Jt_dot(R0t,omega_target)*Jt_inv;


Bcim_dot = [zeros(12,6) zeros(12,6); -Jt_inv_dot*J0-Jt_inv*J0_dot -Jt_inv_dot*Je-Jt_inv*Je_dot];
            
H_star  = Bcim'*H_bar*Bcim;
C1_star = Bcim'*H_bar*Bcim_dot;
C2_star = Bcim'*vertcat(C0,Cq,Ct);
C_star  = C1_star*vertcat(omega0_sat,r0dot_sat,qdot_sat) + C2_star; 

% Rotation matrix to transfer satelite's quantities to the inertial fraùe
% with SPART: omega in body frame, rdot in inertial frame
P0_spart = blkdiag(R0s,eye(3));
P0_dot_spart = blkdiag(P0_spart(1:3,1:3)*SkewSym(omega0_sat),zeros(3));

% Rotation matrix to transfer target's quantities to the inertial frame
% nothing other dans eye(6) because we are already using the dynamics
% expressed in the inertial frame (H_ine, C_ine)
Pt =  blkdiag(eye(3), eye(3));
Pt_dot = blkdiag(Pt(1:3,1:3)*SkewSym(omega_target),zeros(3));

nq = 6; % Number of wheels+joints
Hq_star = H_star(1:nq,1:nq);
Hq0_star = H_star(1:nq,nq+1:end);
H0_star = H_star(nq+1:end,nq+1:end);

Cq_star = C_star(1:nq);
C0_star = C_star(nq+1:end);

% get derivatives for the satelite's inertia matrix (from SPART)
Hdot = robot_pre.Hdot();
H0_dot = Hdot.H0_dotf(R0s,omega0_sat,r0_sat,q_sat,r0dot_sat,qdot_sat);
H0q_dot = Hdot.H0m_dotf(R0s,omega0_sat,r0_sat,q_sat,r0dot_sat,qdot_sat);

% get derivatives for the target's inertia matrix (from free_body_euler)
Htdot = Ht_dot(R0t,omega_target);

% matrices W_0 and W_{q0} in Overleaf
W0 = P0_spart*H0-Pt*Ht*Jt_inv*J0;
W0_dot = P0_dot_spart*H0 + P0_spart*H0_dot - Pt_dot*Ht*Jt_inv*J0 - Pt*Htdot*Jt_inv*J0 - Pt*Ht*Jt_inv_dot*J0 - Pt*Ht*Jt_inv*J0_dot;


Wq0 = P0_spart*H0q-Pt*Ht*Jt_inv*Je;
Wq0_dot = P0_dot_spart*H0q + P0_spart*H0q_dot - Pt_dot*Ht*Jt_inv*Je - Pt*Htdot*Jt_inv*Je - Pt*Ht*Jt_inv_dot*Je -Pt*Ht*Jt_inv*Je_dot;

%simplify what's to come, inverse only once
W0_inv = inv(W0);
W0i_Wq0 = W0_inv*Wq0; 

Qc = [W0_dot Wq0_dot]*vertcat(omega0_sat,r0dot_sat,qdot_sat); 
Bq0 = vertcat(-W0i_Wq0, eye(6));
Bh = vertcat(W0_inv,zeros(6));
Qcbar = vertcat(-W0*Qc,zeros(6,1));

% momentum expressed in the inertial frame
hcapt = Pt*blkdiag(target.I,target.m*eye(3))*vertcat(omega_target, rdot_target);

% computed final matrices
Hdiamond = full(W0i_Wq0'*H0_star*W0i_Wq0 -Hq0_star*W0i_Wq0-W0i_Wq0*Hq0_star'+Hq_star);
Cdiamond = full(Cq_star - W0i_Wq0'*C0_star + W0i_Wq0*W0_inv*[W0_dot Wq0_dot]*(Bq0*vertcat(qdot_sat)+Bh*hcapt));

% final dynamics (\ddot{q} in Overleaf)
ddq = Hdiamond\tau - Hdiamond\Cdiamond;

% equation (13) in Overleaf
ddq0_ddq = Bq0*ddq + Qcbar;
ddq0 = ddq0_ddq(1:6);

% integrate to obtain \dot{q}_0 and \dot{q}, starting from the measured
% quantities
dq0_dq = vertcat(omega0_sat,r0dot_sat,qdot_sat) + 0.1*vertcat(ddq0,ddq);
q0dot = dq0_dq(1:6);
qdot = dq0_dq(7:end);

% equation (5) in Overleaf
ddq0_ddq_ddqt = Bcim*ddq0_ddq + Bcim_dot*dq0_dq;
ddqt = ddq0_ddq_ddqt(13:18);

%integrate to obtain \dot{q}_t
dqt = xtarget(4:end) + 0.1*ddqt;

% gather all results and output to function
ddX = ddq0_ddq_ddqt;
dX = vertcat(dq0_dq,dqt);
end