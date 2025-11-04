function [model]= model_dyn(robot)

%% Eval Mass Matrix

%--- System_dimension ---%
dim_0 = 6;
dim_m = robot.n_q;
dim_f = robot.flexs.nb_mode; 

%--- Initial conditions ---%
R00=eye(3);
r00=zeros(3,1);
qm0=zeros(dim_m,1);
u00=zeros(6,1);
um0=zeros(dim_m,1);


%--- Kinematics ---%
%Kinematics
[~,RL,rJ,rL,e,g]    =   Kinematics(R00,r00,qm0,robot);

%Diferential Kinematics
[Bij,Bi0,P0,pm]     =   DiffKinematics(R00,r00,rL,e,g,robot);
%Velocities
[t0,tL]             =   Velocities(Bij,Bi0,P0,pm,u00,um0,robot);
%Jacobian of the flexible link
[Jf0, Jfm]          =   Jacob_flex(rJ,r00,rL,P0,pm,robot);
[Jf0dot, Jfmdot]    =   Jacobdot_flex(rJ,r00,t0,rL,tL,P0,pm,Jf0, Jfm,u00,um0,robot);
[ J0, Jm, J0d0, Jmd0, J0dm, Jmdm, J0dot, Jmdot, Sx0, Sqm]=Jacob_full(R00,r00,rL,t0,tL,P0,pm,robot);

%--- Dynamics ---%
%Inertias in inertial frames
[I0,Im]             =   I_I(R00,RL,robot);
%Mass Composite Body matrix
[M0_tilde,Mm_tilde, M0, Mm] =   MCB_full(I0,Im,Bij,Bi0,robot);
%Generalized Inertia matrix
[H0,H0m,Hm]         =   GIM(M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot);
[Hf0, Hfm, Hf]      =   GIM_flex(RL,Jf0, Jfm,robot);
%Generalized Convective Inertia matrix
 [C0, C0m, Cm0, Cm] = CIM_rig(tL, J0, Jm, J0d0, Jmd0, J0dm, Jmdm, J0dot, Jmdot, Sx0, Sqm, Mm, u00, um0,robot);
[Cf0, Cfm, Cf]      =   CIM_flex(RL,Jf0dot, Jfmdot,robot);
%Stifness matrix of flexible links
[Kf]                =   Stif_flex(robot);



%% Matrix of the differential second order equation of the system:
% Mm(x).ddot(x)+Dm(x).dot(x)+Km(x).x= w_ext
% with dot(x)=[u0;dot(qm);dot(eta)] and eta the vector of the modal DOF

model.Mm =...
    [H0 H0m Hf0' ;...
    H0m' Hm Hfm';...
    Hf0 Hfm Hf];
model.Dm =...
    [C0 C0m Cf0' ;...
    Cm0 Cm Cfm';...
    Cf0 Cfm Cf];
model.Km=zeros(dim_0+dim_m+dim_f,dim_0+dim_m+dim_f);
model.Km(dim_0+dim_m+1:end,dim_0+dim_m+1:end)=Kf;

end
