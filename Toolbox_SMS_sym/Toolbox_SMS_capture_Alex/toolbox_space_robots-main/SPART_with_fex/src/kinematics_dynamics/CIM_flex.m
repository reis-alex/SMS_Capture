function [Cf0, Cfm, Cf] = CIM_flex(RL,Jf0dot, Jfmdot,robot)
% Computes the flex components of the Generalized Convective Inertia Matrix C of the multibody system.
%
%
% [Cf0, Cfm, Cf] = CIM_flex(RL,Jf0dot, Jfmdot,robot)
%
% :parameters:
%   * RL -- Links CCS 3x3 rotation matrices with respect to the inertial CCS -- as a [3x3xn] matrix.
%   * Jf0dot -- Base-link connection points Jacobian time-derivative -- as a [6x6xn_f] matrix.
%   * Jfmdot -- Manipulator connection points Jacobian time-derivative -- as a [6xn_qxn_f] matrix.
%   * robot -- Robot model (see :doc:`/Tutorial_Robot`).
%
% :return:
%   * Cf0 -- Flexible DOF -- Base-link -- coupling convective inertia matrix -- as a [n_mx6] matrix.
%   * Cfm -- Flexible DOF -- manipulator coupling  convective inertia matrix -- as a [n_mxn_q] matrix.
%   * Cf -- Flexible DOF convective inertia matrix -- as a [n_mxn_m] matrix.
%
% See also: :func:`src.kinematics_dynamics.Jacob_flex`.

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

%--- Number of links and Joints ---%
n_q=robot.n_q;

%--- Number of links  ---%
n_f=robot.n_links_flex;

%--- Number of DOF  ---%
n_m=robot.flexs.nb_mode;


%--- C matrix ---%
%Pre-allocate Cf0
Cf0=zeros(n_m,6,'like',Jf0dot);

%Pre-allocate Cfm
Cfm=zeros(n_m,n_q,'like',Jfmdot);

%Cf
Cf=zeros(n_m,n_m,'like',Jf0dot);


% loop on the flexible links
for ind_l=1:n_f
    %Indice of the flexible DOF
    if ind_l==n_f
        ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.nb_mode;
    else
        ind_DOF=robot.flexs.link_mode_id(ind_l):robot.flexs.link_mode_id(ind_l+1)-1;
    end
    %Indice of the flexible links
    ind_link=robot.flexs.link_id(ind_l);
    
    %rotation matrix from Inertia CCS to the link CCS
    R_CCS_link=[RL(:,:,ind_link)' zeros(3,3,'like',RL); zeros(3,3,'like',RL) RL(:,:,ind_link)'];
    
    Cf0(ind_DOF,:)=robot.flexs.Li(ind_DOF,[4:6 1:3])*R_CCS_link*Jf0dot(:,:,ind_l);
    Cfm(ind_DOF,:)=robot.flexs.Li(ind_DOF,[4:6 1:3])*R_CCS_link*Jfmdot(:,:,ind_l);
    Cf(ind_DOF,ind_DOF)= 2*diag(robot.flexs.pulse(1,ind_DOF).*robot.flexs.damp(1,ind_DOF));
    
end