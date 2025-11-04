function [C0, C0m, Cm0, Cm] =...
    CIM_rig_M(tL, J0, Jm, J0d0, Jmd0, J0dm, Jmdm, J0dot, Jmdot, Sx0, Sqm, Mm, u0, um,robot)
% Computes the Generalized Convective Inertia Matrix C of the multibody system.
%
% :parameters:
%   * tL -- Manipulator twist [\omega,rdot], projected in the inertial CCS -- as a [6xn] matrix.
%   * J0 -- Base-link geometric Jacobian -- as a [6x6xn_l] matrix.
%   * Jm -- Manipulator geometric Jacobian -- as a [6xn_qxn_l] matrix.
%   * J0d0 -- Base-link 2nd order geometric Jacobian -- as a [6x6xn_lx6] matrix.
%   * Jmd0 -- Manipulator 2nd order geometric Jacobian -- as a [6xn_qxn_lx6] matrix.
%   * J0dm -- Base-link 2nd order geometric Jacobian -- as a [6x6xn_lxn_q] matrix.
%   * Jmdm -- Manipulator 2nd order geometric Jacobian -- as a [6xn_qxn_lxn_q] matrix.
%   * J0dot -- Base-link Jacobian time-derivative -- as a [6x6xn_l] matrix.
%   * Jmdot -- Manipulator Jacobian time-derivative -- as a [6xn_qxn_l] matrix.
%   * Sx0 -- Omega0 matrix derivative -- as a [6x6x3] matrix.
%   * Sqm -- Omega matrix derivative -- as a [6x6xn_qxn_l] matrix.
%   * Mm -- Mass and inertia matrices of links projected in the inertial CCS -- as a [3x3xn_l] matrix.
%   * u0 -- Base-link velocities [\omega,rdot]. The angular velocity is projected in the body-fixed CCS,
%           while the linear velocity is projected in the inertial CCS -- [6x1].
%   * um -- Joint velocities -- [n_qx1].
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return: 
%   * C0 -> Base-link convective inertia matrix -- as a [6x6] matrix.
%   * C0m -> Base-link - manipulator coupling convective inertia matrix -- as a [6xn_q] matrix.
%   * Cm0 -> Manipulator - base-link coupling convective inertia matrix -- as a [n_qx6] matrix.
%   * Cm -> Manipulator convective inertia matrix -- as a [n_qxn_q] matrix.

%{  
    LICENSE

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}

%=== CODE ===%

%Pre-allocate Omega
Omega=zeros(6,6,robot.n_links_joints,'like',tL);

%Compute Omega
for i=1:robot.n_links_joints
    Omega(1:6,1:6,i)=[SkewSym(tL(1:3,i)), zeros(3,3);
        zeros(3,3), zeros(3,3)];
end


%Pre-allocate
Idot=zeros(6,6,robot.n_links_joints,'like',tL);

%Compute Idot
for i=1:robot.n_links_joints
    Idot(1:6,1:6,i)=[Omega(1:3,1:3,i)*Mm(1:3,1:3,i), zeros(3,3); zeros(3,6)];
    Idot(1:6,1:6,i)=Idot(1:6,1:6,i)+Idot(1:6,1:6,i)';
end

%--- C Matrix ---%
%Pre-Allocate
C0=zeros(6,6,'like',tL);
Cm=zeros(robot.n_q,robot.n_q,'like',tL);
C0m=zeros(6,robot.n_q,'like',tL);
Cm0=zeros(robot.n_q,6,'like',tL);


for ind_l=1:robot.n_links_joints
    if robot.joints(ind_l).type == 0
        C0=C0+...
            J0(:,:,ind_l).'*Mm(:,:,ind_l)*J0dot(:,:,ind_l)+...
            J0(:,:,ind_l).'*Idot(:,:,ind_l)*J0(:,:,ind_l);
        for ind_0=1:3
            C0(:,ind_0)=C0(:,ind_0)+...
                J0d0(:,:,ind_l,ind_0)'*Mm(:,:,ind_l)*tL(:,ind_l);
        end

    end
    
    if robot.joints(ind_l).type ~= 0

        C0m=C0m+...
            J0(:,:,ind_l).'*Mm(:,:,ind_l)*Jmdot(:,:,ind_l)+...
            J0(:,:,ind_l).'*Idot(:,:,ind_l)*Jm(:,:,ind_l);
        for ind_q=1:robot.n_q
            C0m(:,ind_q)=C0m(:,ind_q)+...
                J0dm(:,:,ind_l,ind_q)'*Mm(:,:,ind_l)*tL(:,ind_l);
        end
        Cm0=Cm0+...
            Jm(:,:,ind_l).'*Mm(:,:,ind_l)*J0dot(:,:,ind_l)+...
            Jm(:,:,ind_l).'*Idot(:,:,ind_l)*J0(:,:,ind_l);
        for ind_0=1:3
            Cm0(:,ind_0)=Cm0(:,ind_0)+...
                Jmd0(:,:,ind_l,ind_0)'*Mm(:,:,ind_l)*tL(:,ind_l);
        end
    
    end
    if (robot.joints(ind_l).type~=0 )

        Cm=Cm+...
            Jm(:,:,ind_l).'*Mm(:,:,ind_l)*Jmdot(:,:,ind_l)+...
            Jm(:,:,ind_l).'*Idot(:,:,ind_l)*Jm(:,:,ind_l);
        for ind_q=1:robot.n_q
        if (robot.joints(ind_l).type~=0 && robot.joints(ind_q).type~=0) && (robot.con.branch(i,ind_q)==1 ...
                || robot.con.branch(ind_q,i)==1)

            Cm(:,ind_q)=Cm(:,ind_q)+...
                Jmdm(:,:,ind_l,ind_q)'*Mm(:,:,ind_l)*tL(:,ind_l);
        end
        end
    end
end

C0d=zeros(6,6,'like',u0);
C0md=zeros(6,robot.n_q,'like',u0);
Cm0d=zeros(robot.n_q,6,'like',u0);
Cmd=zeros(robot.n_q,robot.n_q,'like',u0);

tLdx0=zeros(6,robot.n_links_joints,6,'like',tL);


for ind_0=1:6
    for ind_l=1:robot.n_links_joints
        tLdx0(:,ind_l,ind_0)=...
            J0d0(:,:,ind_l,ind_0)*u0+Jmd0(:,:,ind_l,ind_0)*um;
        C0d(ind_0,:)=C0d(ind_0,:)+...
            tLdx0(:,ind_l,ind_0)'*Mm(:,:,ind_l)*J0(:,:,ind_l)+...
            tL(:,ind_l)'*Mm(:,:,ind_l)*J0d0(:,:,ind_l,ind_0);
        C0md(ind_0,:)=C0md(ind_0,:)+...
            tLdx0(:,ind_l,ind_0)'*Mm(:,:,ind_l)*Jm(:,:,ind_l)+...
            tL(:,ind_l)'*Mm(:,:,ind_l)*Jmd0(:,:,ind_l,ind_0);
        if ind_0<4
            dMm= Sx0(:,:,ind_0)*Mm(:,:,ind_l)+Mm(:,:,ind_l)*Sx0(:,:,ind_0)';
            C0d(ind_0,:)=C0d(ind_0,:)+...
                tL(:,ind_l)'*dMm*J0(:,:,ind_l);
            C0md(ind_0,:)=C0md(ind_0,:)+...
                tL(:,ind_l)'*dMm*Jm(:,:,ind_l);
        end
    end
end

tLdxm=zeros(6,robot.n_links_joints,robot.n_q,'like',tL);


for ind_m=1:robot.n_q
    for ind_l=1:robot.n_links_joints
        tLdxm(:,ind_l,ind_m)=...
            J0dm(:,:,ind_l,ind_m)*u0+Jmdm(:,:,ind_l,ind_m)*um;
        dMm= Sqm(:,:,ind_m,ind_l)*Mm(:,:,ind_l)+Mm(:,:,ind_l)*Sqm(:,:,ind_m,ind_l)';
        Cm0d(ind_m,:)=Cm0d(ind_m,:)+...
            tLdxm(:,ind_l,ind_m)'*Mm(:,:,ind_l)*J0(:,:,ind_l)+...
            tL(:,ind_l)'*Mm(:,:,ind_l)*J0dm(:,:,ind_l,ind_m);
        Cmd(ind_m,:)=Cmd(ind_m,:)+...
            tLdxm(:,ind_l,ind_m)'*Mm(:,:,ind_l)*Jm(:,:,ind_l)+...
            tL(:,ind_l)'*...
            (dMm*Jm(:,:,ind_l)+Mm(:,:,ind_l)*Jmdm(:,:,ind_l,ind_m));
    end
end

C0=C0-C0d/2;
C0m=C0m-C0md/2;
Cm0=Cm0-Cm0d/2;
Cm=Cm-Cmd/2;