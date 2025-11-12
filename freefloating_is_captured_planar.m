function [xsat_plus, xtarget_plus] = freefloating_is_captured_planar(dt, xsat, R0s, xtarget, R0t, robot_pre, robot_pre_sim, target, tau, hcapt)
tic
% load state vector for satelite and target
r0_sat          = xsat(satelite.idx.r0);
q_sat           = xsat(satelite.idx.q);
omega0_sat      = xsat(satelite.idx.omega0);
r0dot_sat       = xsat(satelite.idx.r0dot);
qdot_sat        = xsat(satelite.idx.qdot);

pos_target   = xtarget(target.idx.r0,1);
omega_target = xtarget(target.idx.omega0,1);
rdot_target  = xtarget(target.idx.r0dot,1);

% Rotation matrix to transfer satelite's quantities to the inertial fraùe
% with SPART: omega in body frame, rdot in inertial frame
P0_spart = blkdiag(R0s,eye(3));
Pt =  blkdiag(R0t, eye(3));


% get matrices from the satelites dynamics (from SPART)
H   = full(satelite.dynamics.H(R0s,xsat));
Hq  = H(7:12,7:12);
H0q = H(1:6,7:12);
H0  = H(1:6,1:6);

% get matrices from the target dynamics (from SPART)
Ht      = full(target.dynamics.H(R0t,xtarget));
H_bar   = full(blkdiag(H,Ht));

% end effector jacobian (from SPART)
EE_idx = find(strcmp({robot_pre_sim.links.name},'Link_EE'));
Jee = satelite.Jacob(EE_idx);
Je  = full(Jee.Jmf(R0s,xsat));
J0  = full(Jee.J0f(R0s,xsat));

% target jacobian

% Jt      = target.Jacob(1); 
% Jt      = full(Jt.J0f(R0t,pos_target,0));
Jt = blkdiag(R0t,eye(3));
Jt_inv  = blkdiag(R0t',eye(3)); % after discussion with Mathieu (Jt computed wrt target's CoM, numerically better)
                                % add SkewSym(rt-rp), bras de levier

%% initialize capture


% moment = [[P0_spart*H0 P0_spart*H0q]*xsat(13:end) + Pt*Ht*xtarget(7:end)];

% From here and on, all matrices are computed with the constrained
% velocities

C = full(robot_pre.dynamics.C(R0s,r0_sat,q_sat,omega0_sat,r0dot_sat,qdot_sat))*[omega0_sat;r0dot_sat;qdot_sat];
C0 = C(1:6);
Cq = C(7:end);

Ct      = full(full(target.dynamics.C(R0t,pos_target,0,omega_target,rdot_target,0))*[omega_target; rdot_target]);

% compute Je_dot and J0_dot from robot_pre
Jee_dot = robot_pre.Jacobdot(EE_idx);
Je_dot = Jee_dot.Jmdotf(R0s,omega0_sat,r0_sat,q_sat,r0dot_sat,qdot_sat);
J0_dot = Jee_dot.J0dotf(R0s,omega0_sat,r0_sat,q_sat,r0dot_sat,qdot_sat);

% calc d(Jt_inv)/dt
Jt_dot = target.Jacobdot(1);
Jt_inv_dot = full(-Jt_inv*Jt_dot.J0dotf(R0t,omega_target,pos_target,0,rdot_target,0)*Jt_inv);


% Computed Bcim and Bcim_dot
Bcim = full([eye(6,6)     zeros(6,3);
            zeros(3,6)    eye(3);
            -Jt_inv*J0 -Jt_inv*Je]);
Bcim_dot = full([zeros(9,6) zeros(9,3); -Jt_inv_dot*J0-Jt_inv*J0_dot -Jt_inv_dot*Je-Jt_inv*Je_dot]);
       
% Compute Hstar and Cstar
H_star  = full(Bcim'*H_bar*Bcim);
C1_star = full(Bcim'*H_bar*Bcim_dot);
C2_star = Bcim'*vertcat(C0,Cq,Ct);
C_star  = full(C1_star*vertcat(omega0_sat,r0dot_sat,qdot_sat) + C2_star); 

% get derivatives for the satelite's inertia matrix (from SPART)

Hdot    = robot_pre.dynamics.Hdot(R0s,omega0_sat,r0_sat,q_sat,r0dot_sat,qdot_sat);
H0_dot  = Hdot(1:6,1:6);
H0q_dot = Hdot(1:6,7:end);
% get derivatives for the target's inertia matrix (from SPART)

Htdot = full(target.dynamics.Hdot(R0t,omega_target,pos_target,0,rdot_target,0));

% compute derivatives of P0 and Pt
P0_dot_spart = blkdiag(R0s*SkewSym(omega0_sat),zeros(3));
Pt_dot = blkdiag(R0t*SkewSym(omega_target),zeros(3));

% matrices W_0 and W_{q0} in Overleaf
W0 = P0_spart*H0-Pt*Ht*Jt_inv*J0;
W0_dot = P0_dot_spart*H0 + P0_spart*H0_dot - Pt_dot*Ht*Jt_inv*J0 - Pt*Htdot*Jt_inv*J0 - Pt*Ht*Jt_inv_dot*J0 - Pt*Ht*Jt_inv*J0_dot;

Wq0 = P0_spart*H0q-Pt*Ht*Jt_inv*Je;
Wq0_dot = P0_dot_spart*H0q + P0_spart*H0q_dot - Pt_dot*Ht*Jt_inv*Je - Pt*Htdot*Jt_inv*Je - Pt*Ht*Jt_inv_dot*Je -Pt*Ht*Jt_inv*Je_dot;

%simplify what's to come, inverse only once
W0_inv = inv(W0);
W0i_Wq0 = full(W0\Wq0); 
W0_inv_hcapt = full(W0\hcapt);

Qc = [W0_dot Wq0_dot]*vertcat(omega0_sat,r0dot_sat,qdot_sat);
Bq0 = vertcat(-W0i_Wq0, eye(3,3));
Bh = vertcat(W0_inv_hcapt, zeros(3,1));
Qcbar = vertcat(- (W0 \ Qc), zeros(3,1));

% computed final matrices
Hdiamond = Bq0'*H_star*Bq0;
Cdiamond = Bq0'*(C_star+full(Qcbar));

xsat_plus = zeros(18,1);
xtarget_plus = zeros(12,1);

% final dynamics (\ddot{q} in Overleaf)
ddq = Hdiamond \ (tau - Cdiamond);

% integrate ddq
qdot_plus = qdot_sat + dt*full(ddq);
qdot0_plus = - (W0 \ (Wq0 * qdot_plus)) + (W0 \ hcapt);
xsat_plus(16:18) = qdot_plus;
xsat_plus(10:15) = qdot0_plus;

qdot_full_post = vertcat(xsat_plus(10:15), xsat_plus(16:end));
xtarget_plus(7:12) = (-Jt_inv * [J0, Je]) * qdot_full_post; 
xsat_plus(1:9) = xsat(1:9) + dt * xsat_plus(10:end);   % use xsat_plus velocities (13:24)
xtarget_plus(1:6) = xtarget(1:6) + dt * xtarget_plus(7:12);
% compute dq0 from Equation (12 in Overleaf)
% xsat_plus(13:18) = -W0i_Wq0*xsat_plus(19:end) + Bh(1:6,:)*hcapt;
xsat_plus(13:18) = -W0i_Wq0*xsat(16:end) + W0_inv*hcapt;

% compute dqt equation (4) in Overleaf
% xtarget_plus(7:end) = Bcim(13:18,:)*xsat_plus(13:end);
% xtarget_plus(7:end) = [-Jt_inv*J0 -Jt_inv*Je]*xsat(13:end);
% 
% % integrate velocities to obtain positions, output them
% xsat_plus(1:12) = xsat(1:12) + dt*xsat(13:end);
% xtarget_plus(1:6) = xtarget(1:6) + dt*xtarget(7:end);

end