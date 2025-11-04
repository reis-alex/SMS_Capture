function [C0, C0m, Cm0, Cm] = CIM_rigid_bis(t0,tL,I0,Im,M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,u0,um,robot)
% Computes the Generalized Convective Inertia Matrix C of the multibody system.
%
% :parameters: 
%   * t0 -- Base-link twist [\omega,rdot], projected in the inertial CCS -- as a [6x1] matrix.
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * I0 -- Base-link inertia matrix, projected in the inertial CCS -- as a [3x3] matrix.
%   * Im -- Links inertia matrices, projected in the inertial CCS -- as a [3x3xn] matrix.
%   * M0_tilde -- Base-link mass composite body matrix -- as a [6x6] matrix .
%   * Mm_tilde -- Manipulator mass composite body matrix -- as a [6x6xn] matrix.
%   * Bij -- Twist-propagation matrix (for manipulator i>0 and j>0) -- as a [6x6xn] matrix.
%   * Bi0 -- Twist-propagation matrix (for i>0 and j=0) -- as a [6x6xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * C0 -> Base-link convective inertia matrix -- as a [6x6] matrix.
%   * C0m -> Base-link - manipulator coupling convective inertia matrix -- as a [6xn_q] matrix.
%   * Cm0 -> Manipulator - base-link coupling convective inertia matrix -- as a [n_qx6] matrix.
%   * Cm -> Manipulator convective inertia matrix -- as a [n_qxn_q] matrix.
%
% To obtain the full convective inertia matrix C:
%
% .. code-block:: matlab
%   
%   %Compute the Convective Inertia Matrix C
%   [C0, C0m, Cm0, Cm] = CIM(t0,tL,I0,Im,M0_tilde,Mm_tilde,Bij,Bi0,P0,pm,robot)
%   C=[C0,C0m;Cm0,Cm];
%
% See also: :func:`src.kinematics_dynamics.GIM`.


%=== CODE ===%

%--- Number of links and Joints ---%
n_q=robot.n_q;
n=robot.n_links_joints;

%--- Omega ---%
%Base-link Omega
Omega0=[SkewSym(t0(1:3)), zeros(3,3);
    zeros(3,3), zeros(3,3)];

%Pre-allocate Omega
Omega=zeros(6,6,n,'like',t0);

%Compute Omega
for i=1:n
    Omega(1:6,1:6,i)=[SkewSym(tL(1:3,i)), zeros(3,3);
        zeros(3,3), SkewSym(tL(1:3,i))];
end

%--- Mdot ---%
%Base-link Mdot
Mdot0=[Omega0(1:3,1:3)*I0, zeros(3,3); zeros(3,3), zeros(3,3)];

%Pre-allocate
Mdot=zeros(6,6,n,'like',t0); 

%Compute Mdot
for i=1:n
    Mdot(1:6,1:6,i)=[Omega(1:3,1:3,i)*Im(1:3,1:3,i), zeros(3,3); zeros(3,3), zeros(3,3)];
end

%--- Mdot tilde ---%
%Pre-Allocate
Mdot_tilde=zeros(6,6,n,'like',t0);

%Backwards recursion
for i=n:-1:1
    %Initialize
    Mdot_tilde(1:6,1:6,i)=Mdot(1:6,1:6,i);
    %Add children contributions
    child=find(robot.con.child(:,i))';
    for j=1:length(child)
        Mdot_tilde(1:6,1:6,i)=Mdot_tilde(1:6,1:6,i)+Mdot_tilde(1:6,1:6,child(j));
    end
end
%Base-link
Mdot0_tilde = Mdot0;
%Add children contributions
child=find(robot.con.child_base)';
for j=1:length(child)
    Mdot0_tilde = Mdot0_tilde+Mdot_tilde(1:6,1:6,child(j));
end


%--- Bdot ---%

%Pre-allocate Bdotij
Bdotij=zeros(6,6,n,n,'like',t0);

%Compute Bdotij
for j=1:n
    for i=1:n
        if robot.con.branch(i,j)==1
            %Links are in the same branch
            Bdotij(1:6,1:6,i,j)=[zeros(3,3), zeros(3,3); SkewSym(tL(4:6,j)-tL(4:6,i)), zeros(3,3)];
        else
            %Links are not in the same branch
            Bdotij(1:6,1:6,i,j)=zeros(6,6);
        end
    end
end


%--- Hij tilde ---%
%Pre-allocate Hij_tilde
Hij_tilde=zeros(6,6,n,n,'like',t0);

%Hij_tilde
for i=n:-1:1
    for j=n:-1:1
        Hij_tilde(1:6,1:6,i,j)=Mm_tilde(1:6,1:6,i)*Bdotij(1:6,1:6,i,j);
        %Add children contributions
        child=find(robot.con.child(:,i))';
        for k=1:length(child)
            Hij_tilde(1:6,1:6,i,j)=Hij_tilde(1:6,1:6,i,j)+Bij(1:6,1:6,child(k),i)'*Hij_tilde(1:6,1:6,child(k),i);
        end
    end
end

%Pre-allocate Hi0_tilde and H0j_tilde
Hi0_tilde=zeros(6,6,n,'like',t0);

%Hi0_tilde
for i=n:-1:1
    Bdot=[zeros(3,3), zeros(3,3); SkewSym(t0(4:6)-tL(4:6,i)), zeros(3,3)];
    Hi0_tilde(1:6,1:6,i)=Mm_tilde(1:6,1:6,i)*Bdot;
    %Add children contributions
    child=find(robot.con.child(:,i))';
    for k=1:length(child)
        Hi0_tilde(1:6,1:6,i)=Hi0_tilde(1:6,1:6,i)+Bij(1:6,1:6,child(k),i)'*Hij_tilde(1:6,1:6,child(k),i);
    end
end

% %--- C Matrix ---%
% %Pre-allocate
% Cm=zeros(n_q,n_q,'like',t0);
% C0m=zeros(6,n_q,'like',t0);
% Cm0=zeros(n_q,6,'like',t0);
% 
% %Cm Matrix
% for j=1:n
%     for i=1:n
%         %Joints must not be fixed and links on the same branch
%         if (robot.joints(i).type~=0 && robot.joints(j).type~=0) && (robot.con.branch(i,j)==1 || robot.con.branch(j,i)==1)
%             %Compute Cm matrix
%             if i<=j
%                 %Add children contributions
%                 child_con=zeros(6,6);
%                 child=find(robot.con.child(:,j))';
%                 for k=1:length(child)
%                     child_con=child_con+Bij(1:6,1:6,child(k),i)'*Hij_tilde(1:6,1:6,child(k),j);
%                 end
%                 Cm(robot.joints(i).q_id,robot.joints(j).q_id)=pm(1:6,i)'*(Bij(1:6,1:6,j,i)'*Mm_tilde(1:6,1:6,j)*Omega(1:6,1:6,j)...
%                     +child_con+Mdot_tilde(1:6,1:6,j))*pm(1:6,j);
%             else
%                 Cm(robot.joints(i).q_id,robot.joints(j).q_id)=pm(1:6,i)'*(Mm_tilde(1:6,1:6,i)*Bij(1:6,1:6,i,j)*Omega(1:6,1:6,j)...
%                     +Hij_tilde(1:6,1:6,i,j)+Mdot_tilde(1:6,1:6,i))*pm(1:6,j);
%             end
%         end
%     end
% end
% %C0 matrix
% %Add children contributions
% child_con=zeros(6,6);
% child=find(robot.con.child_base)';
% for k=1:length(child)
%     child_con=child_con+Bi0(1:6,1:6,child(k))'*Hi0_tilde(1:6,1:6,child(k));
% end
% C0 = P0'*(M0_tilde*Omega0+child_con+Mdot0_tilde)*P0;
% 
% %C0m
% for j=1:n
%     if  robot.joints(j).type~=0
%         if j==n
%          C0m(1:6,robot.joints(j).q_id)=P0'*(Bi0(1:6,1:6,j)'*Mm_tilde(1:6,1:6,j)*Omega(1:6,1:6,j)+Mdot_tilde(1:6,1:6,j))*pm(1:6,j);
%         else
%             %Add children contributions
%             child_con=zeros(6,6);
%             child=find(robot.con.child(:,j))';
%             for k=1:length(child)
%                 child_con=child_con+Bi0(1:6,1:6,child(k))'*Hij_tilde(1:6,1:6,child(k),j);
%             end
%             C0m(1:6,robot.joints(j).q_id)=P0'*(Bi0(1:6,1:6,j)'*Mm_tilde(1:6,1:6,j)*Omega(1:6,1:6,j)+child_con...
%                 +Mdot_tilde(1:6,1:6,j))*pm(1:6,j);
%         end
%     end
% end
% %Cm0
% for i=1:n
%     if  robot.joints(i).type~=0  
%         Cm0(robot.joints(i).q_id,1:6)=pm(1:6,i)'*(Mm_tilde(1:6,1:6,i)*Bi0(1:6,1:6,i)*Omega0+Hi0_tilde(1:6,1:6,i)...
%             +Mdot_tilde(1:6,1:6,i))*P0;
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--Compute dT/dx0 and dT/dq --%

H0dot   = zeros(6,6,'like',t0);
dH0dx0  = zeros(6,6,'like',t0);
dH0mdx0 = zeros(6,n_q,'like',t0);
dH0mdq  = zeros(6,n_q,'like',t0);
dHmdq  = zeros(n_q,n_q,'like',t0);


P0dot = Omega0*P0;

%--Compute H0dot --%

H0dot = P0dot'*M0_tilde*P0 + P0'*M0_tilde*P0dot + P0' * Mdot0_tilde *P0;

%--Compute H0mdot --%
H0mdot = zeros(6,n_q,'like',t0);

for i=1:n
    if robot.joints(i).type~=0
        Bi0dot = [zeros(3,6);SkewSym(t0(4:6)-tL(4:6,i)) zeros(3)];
        
        H0mdot(1:6,robot.joints(i).q_id) = P0dot.' * Bi0(1:6,1:6,i).' *  Mm_tilde(1:6,1:6,i).' * pm(1:6,i) +...
                                           P0.' * Bi0dot.' * Mm_tilde(1:6,1:6,i).' * pm(1:6,i) +...
                                           P0.' * Bi0(1:6,1:6,i).' *  Mdot_tilde(1:6,1:6,i).' * pm(1:6,i);
    end
end

H0mdotT = H0mdot.';

%--Compute Hmdot --%
Hmdot = zeros(n_q,n_q);
for j=1:n
    for i=j:n
        if robot.joints(i).type~=0 && robot.joints(j).type~=0
        Hmdot(robot.joints(i).q_id,robot.joints(j).q_id) = pm(1:6,i).' * Mdot_tilde(1:6,1:6,i).' * Bij(1:6,1:6,i,j).' *pm(1:6,j)+...
                                                           pm(1:6,i).' * Mm_tilde(1:6,1:6,i).' * Bdotij(1:6,1:6,i,j).' *pm(1:6,j); 
            
            Hmdot(robot.joints(j).q_id,robot.joints(i).q_id)=Hmdot(robot.joints(i).q_id,robot.joints(j).q_id);
        end
    end
end

%compute dHm/dx0
dHmdx0 = zeros(6,n_q,'like',t0);

%compute dH0/dx0

dM0_tilde_dx0 = zeros(6,6);
dBi0dx0 = [zeros(3) zeros(3);eye(3) zeros(3)]; 
child=find(robot.con.child_base)';   %Add children contributions
for j=1:length(child)
    dMm_tilde_dx0 = zeros(6,6);       %% Bij ne depends pas de x0...
    dM0_tilde_dx0 = dM0_tilde_dx0 + dBi0dx0.' * Mm_tilde(1:6,1:6,child(j)) * Bi0(1:6,1:6,child(j)) +...
                    Bi0(1:6,1:6,child(j))' * dMm_tilde_dx0 *Bi0(1:6,1:6,child(j)) +...
                    Bi0(1:6,1:6,child(j))' * Mm_tilde(1:6,1:6,child(j)) * dBi0dx0;
end
dP0dx0 = [eye(3) zeros(3);zeros(3,6)]; 
dBi0dx0 = [zeros(3) zeros(3);eye(3) zeros(3)];

 
dH0dx0 = dP0dx0.' * M0_tilde * P0 + P0.' * dM0_tilde_dx0 * P0 + P0.' * M0_tilde * dP0dx0;


%compute dH0m/dx0
dH0mdx0 = zeros(6,n_q);
for i=1:n
    if robot.joints(i).type~=0
        dH0mdx0(1:6,robot.joints(i).q_id) = dP0dx0.' * Bi0(1:6,1:6,i).' * Mm_tilde(1:6,1:6,i).' * pm(1:6,i) + ...
                                           P0.'  * dBi0dx0.' * Mm_tilde(1:6,1:6,i).' * pm(1:6,i) ;
    end
end

%compute dH0/dq
dH0dq = zeros(6,6,'like',t0);

%compute dH0m/dq

dH0mdq = zeros(6,n_q);
dBijdqi = [zeros(3) zeros(3);eye(3) zeros(3)];
dBi0dqi = [zeros(3) zeros(3);-eye(3) zeros(3)];
for i=1:n
    if robot.joints(i).type~=0
        dMm_tilde_dq = zeros(6,6);
        child=find(robot.con.child(:,i))';
        for j=1:length(child)
            dMm_tilde_dq = dMm_tilde_dq + dBijdqi'*Mm_tilde(1:6,1:6,child(j))*Bij(1:6,1:6,child(j),i) ...
                           + Bij(1:6,1:6,child(j),i).' *Mm_tilde(1:6,1:6,child(j)) *dBijdqi;
        end
        
        dH0mdq(1:6,robot.joints(i).q_id) = P0.' * dBi0dqi.' * Mm_tilde(1:6,1:6,i).' * pm(1:6,i) +...
                                          P0.' * Bi0(1:6,1:6,i).' * dMm_tilde_dq.' * pm(1:6,i);
    end
end


%compute dHm/dq
dMm_tilde_dq = zeros(6,6,'like',t0);
for i=1:n
    if robot.joints(i).type~=0
        dMm_tilde_dq = zeros(6,6);
        child=find(robot.con.child(:,i))';
        for j=1:length(child)
            dMm_tilde_dq = dMm_tilde_dq + dBijdqi'*Mm_tilde(1:6,1:6,child(j))*Bij(1:6,1:6,child(j),i) ...
                           + Bij(1:6,1:6,child(j),i).' *Mm_tilde(1:6,1:6,child(j)) *dBijdqi;
        end
    end
end

for j=1:n
    for i=j:n
        if robot.joints(i).type~=0 && robot.joints(j).type~=0
            
            dHmdq(robot.joints(i).q_id,robot.joints(j).q_id)  = pm(1:6,i).' * dMm_tilde_dq.' * Bij(1:6,1:6,i,j).'*pm(1:6,j) +...
                                                                pm(1:6,i).' * Mm_tilde(1:6,1:6,i) * dBijdqi.' *pm(1:6,j) ;
       
            dHmdq(robot.joints(j).q_id,robot.joints(i).q_id)=dHmdq(robot.joints(i).q_id,robot.joints(j).q_id);
        end
    end
end


%--- C Matrix ---%
%Pre-allocate
Cm=zeros(n_q,n_q,'like',t0);
C0=zeros(6,6,'like',t0);
C0m=zeros(6,n_q,'like',t0);
Cm0=zeros(n_q,6,'like',t0);


C0 = H0dot - 0.5 * diag(u0.' * dH0dx0) - 0.5 * diag(um.' * dH0mdx0.');
%Cm = Hmdot - 0.5 * um.' * dHmdq - 0.5 * diag(u0.' * dH0mdq); %og
Cm = Hmdot - 0.5 * diag(um.') * dHmdq - 0.5 * diag(u0.' * dH0mdq);
C0m = H0mdot - 0.5 * (diag(u0.') * dH0mdx0);% -  0.5 * diag(um.' * dHmdx0) ;
Cm0 = H0mdotT - 0.5 * diag(um.') * dH0mdq.';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

