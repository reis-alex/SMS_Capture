function [ J0, Jm, J0d0, Jmd0, J0dm, Jmdm, J0dot, Jmdot, Sx0, Sqm]=...
    Jacob_full(R0,r0,rL,t0,tL,P0,pm,robot)
% Computes Jacobians of links COM.
%
% [ J0, Jm, J0d0, Jmd0, J0dm, Jmdm, J0dot, Jmdot, Sx0, Sqm]=Jacob_full(R0,r0,rL,t0,tL,P0,pm,robot)
%
% :parameters:
%   * r0 -- Position of the base-link center-of-mass with respect to the origin of the inertial frame,
%projected in the inertial CCS -- [3x1].
%   * rL -- Positions of the links, projected in the inertial CCS -- as a [3xn] matrix.
%   * P0 -- Base-link twist-propagation "vector" -- as a [6x6] matrix.
%   * pm -- Manipulator twist-propagation "vector" -- as a [6xn] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return:
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

%Pre-allocate


J0=zeros(6,6,robot.n_links_joints,'like',r0);
Jm=zeros(6,robot.n_q,robot.n_links_joints,'like',r0);

J0dot=zeros(6,6,robot.n_links_joints,'like',r0);
Jmdot=zeros(6,robot.n_q,robot.n_links_joints,'like',r0);

J0d0=zeros(6,6,robot.n_links_joints,6,'like',r0);
Jmd0=zeros(6,robot.n_q,robot.n_links_joints,6,'like',r0);

J0dm=zeros(6,6,robot.n_links_joints,robot.n_q,'like',r0);
Jmdm=zeros(6,robot.n_q,robot.n_links_joints,robot.n_q,'like',r0);

%Jacobian of the links COM

for ind=1:robot.n_links_joints
    [J0(:,:,ind), Jm(:,:,ind)]=Jacob(rL(:,ind),r0,rL,P0,pm,ind,robot);
    [J0dot(:,:,ind), Jmdot(:,:,ind)]=Jacobdot(rL(:,ind),tL(:,ind),r0,t0,rL,tL,P0,pm,ind,robot);
end


%Pre-allocate
Sx0=zeros(6,6,3,'like',r0);
Sq0=zeros(6,6,3,'like',r0);
Sqm=zeros(6,6,robot.n_q,robot.n_links_joints,'like',r0);

for i=1:3
    Sx0(1:3,1:3,i)=SkewSym(R0(:,i));
    Sq0(1:3,1:3,i)=Sx0(1:3,1:3,i);
    Sq0(4:6,4:6,i)=Sx0(1:3,1:3,i);
end

for ind_q=1:robot.n_q
    for ind_j=1:robot.n_links_joints
        Sqm(1:3,1:3,ind_q,ind_j)=SkewSym(Jm(1:3,ind_q,ind_j));
        Sqm(4:6,4:6,ind_q,ind_j)=Sqm(1:3,1:3,ind_q,ind_j);
    end
end

for ind_j=1:robot.n_links_joints
    
    %Jacobian derivative relative to \omega
    for ind_0=1:3
        J0d0(:,:,ind_j,ind_0)=[eye(3),zeros(3,3);SkewSym(r0-rL(:,ind_j)),eye(3)]*Sx0(:,:,ind_0)*P0+...
            [zeros(3,3),zeros(3,3);SkewSym(-J0(4:6,ind_0,ind_j)),zeros(3,3)]*P0;
        for ind_q=1:ind_j
            %If joint is not fixed
            if (robot.joints(ind_q).type~=0) && (robot.con.branch(ind_j,ind_q)==1)
                Jmd0(:,robot.joints(ind_q).q_id,ind_j,ind_0)=...
                    [eye(3),zeros(3,3);SkewSym(rL(1:3,ind_q)-rL(1:3,ind_j)),eye(3)]*Sq0(:,:,ind_0)*pm(1:6,ind_q)+...
                    [zeros(3,6);SkewSym(J0(4:6,ind_0,ind_q)-J0(4:6,ind_0,ind_j)),eye(3)]*pm(1:6,ind_q);
            end
        end
    end
    
    %Jacobian derivative relative to \qm
    for ind_q_d=1:ind_j
        if (robot.joints(ind_q_d).type~=0) && (robot.con.branch(ind_j,ind_q_d)==1)
            q_d_id=robot.joints(ind_q_d).q_id;
            J0dm(:,:,ind_j,q_d_id)=...
                [zeros(3,6);SkewSym(-Jm(4:6,q_d_id,ind_j)),zeros(3,3)]*P0;
            for ind_q=1:ind_q_d
                if (robot.joints(ind_q).type~=0) && (robot.con.branch(ind_j,ind_q)==1)
                    q_id=robot.joints(ind_q).q_id;
                    Jmdm(:,q_id,ind_j,q_d_id)=...
                    [eye(3),zeros(3,3);SkewSym(rL(1:3,ind_q)-rL(1:3,ind_j)),eye(3)]*Sqm(:,:,q_d_id,ind_j)*pm(1:6,ind_q)+...
                        [zeros(3,6);SkewSym(Jm(4:6,q_d_id,q_id)-Jm(4:6,q_d_id,ind_j)),zeros(3,3)]*pm(1:6,ind_q);
                end
            end
            
        end
    end
end