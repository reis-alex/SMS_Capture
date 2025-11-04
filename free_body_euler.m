classdef free_body_euler
    properties
        m
        I
        angles
        b_omega
        b_vel
        i_pos
        i_vel
        b_tau
        R
    end

    methods
        function obj = free_body_euler(params)
            obj.I = params.I;
            obj.m = params.m;
%             obj.angles = casadi.SX.sym('angles',3,1);
            obj.R = casadi.SX.sym('rotation_matrix',3,3);
            obj.b_omega = casadi.SX.sym('body_omega',3,1);
            obj.b_vel = casadi.SX.sym('body_velocity',3,1);
            obj.i_vel = casadi.SX.sym('inertial_velocity',3,1);
            obj.i_pos = casadi.SX.sym('inertial_position',3,1);
            obj.b_tau = casadi.SX.sym('body_wrench',6,1);
        end

        function out = SkewSym(obj,a)
            out = [ 0   -a(3)  a(2);
                   a(3)   0   -a(1);
                  -a(2)  a(1)   0 ];
        end

        function out = T_IB(obj)
            phi = obj.angles(1); theta = obj.angles(2); psi = obj.angles(3);
            cphi = cos(phi); sphi = sin(phi);
            cth  = cos(theta); sth = sin(theta);
            cpsi = cos(psi); spsi = sin(psi);

            Rz = [ cpsi -spsi 0;
                spsi  cpsi 0;
                0     0    1];
            Ry = [ cth 0 sth;
                0   1  0;
                -sth 0 cth];
            Rx = [1 0 0;
                0 cphi -sphi;
                0 sphi  cphi];
            Rot = Rz * Ry * Rx; % body -> inertial
            out = blkdiag(Rot, eye(3)); % velocity already in inertial frame
        end

        function out = T_IBf(obj)
            TIB = obj.T_IB();
            out = casadi.Function('T_body_to_inertial',{obj.angles},{TIB},{'Angles'},{'T_IB'});
        end

        function out = T_IB_dot(obj)
            Rot = obj.R;
            Rot = Rot(1:3,1:3);
            Rot_dot = Rot*obj.SkewSym(obj.b_omega);
            out = blkdiag(Rot_dot,zeros(3));
        end

        function out = body_H(obj)
            out = blkdiag(obj.I, obj.m*eye(3));
        end

        function out = body_C(obj)
            C = casadi.SX.zeros(6,6);
            C(4:6,4:6) = zeros(3);% obj.m*obj.SkewSym(obj.b_omega);    % Coriolis for linear part
            C(1:3,1:3) = obj.SkewSym(obj.I*obj.b_omega);      % Gyroscopic part %sign?
            out = casadi.Function('C',{obj.b_omega},{C},{'Omega'},{'C'});
        end

        function out = ine_H(obj)
            H_body = obj.body_H();
            T_IB = blkdiag(obj.R,eye(3));
            H_ine = T_IB'*(H_body*T_IB);
            out = casadi.Function('H_inertial',{obj.R},{H_ine},{'Angles'},{'H_inertial'});
        end

        function out = ine_H_dot(obj)
            H_body = obj.body_H();
            T_IB = blkdiag(obj.R,eye(3));
            dT_IB = blkdiag(obj.R*SkewSym(obj.b_omega), zeros(3));

            T_inv = inv(T_IB);
            T_invT = T_inv';

            Hdot = T_invT * ( -dT_IB' * H_body - H_body * dT_IB ) * T_inv;

            out = casadi.Function('Hdot',...
                {obj.R, obj.b_omega}, {Hdot},...
                {'RotationMatrix','Omega'},{'Hdot'});

        end

        function out = ine_C(obj)
            H_body = obj.body_H();
            C_body = obj.body_C();
            T_IB = blkdiag(obj.R,eye(3));
            dT_IB = blkdiag(obj.R*SkewSym(obj.b_omega), zeros(3));

            Cine = T_IB'*C_body(obj.b_omega);
            out = casadi.Function('C_inertial',{obj.R, obj.b_omega},{Cine},{'RotationMatrix','Omega'},{'C_inertial'});
        end

        function out = J_point_body(obj,r_b)
            out = [ eye(3), zeros(3);
                    -obj.SkewSym(r_b), eye(3) ];
        end

        function out = J_point_inertial(obj,r_b)
            AdR = blkdiag(obj.R,eye(3));
            Jb = obj.J_point_body(r_b);
            Ji = AdR * Jb;
            out = casadi.Function('Jine',{obj.R},{Ji},{'RotationMatrix'},{'Jacobian_to_ine'});
        end

        % future: simple inverse of J_point_inertial to be coded

        function out = Jacobdot(obj,r_b)
            Rot = obj.R;
            Rot = Rot(1:3,1:3);
            Rot_dot = Rot*obj.SkewSym(obj.b_omega);
            Jdot = [Rot_dot zeros(3,3); -Rot_dot*obj.SkewSym(obj.b_omega)*obj.SkewSym(r_b) Rot_dot];
            out = casadi.Function('Jacob_dot',{obj.R,obj.b_omega},{Jdot},{'RotationMatrix','omega_body'},{'Jacob_dot'});
        end

        function out = ddX_ine(obj)
            H_ine = obj.ine_H();
            C_ine = obj.ine_C();
            dq = [obj.b_omega; obj.i_vel];
            T_IB = blkdiag(obj.R,eye(3));
            ddX = H_ine(obj.R) \ (T_IB*obj.b_tau - C_ine(obj.R, obj.b_omega) * dq);

            out = casadi.Function('ddX_inertial', ...
                {obj.R, obj.b_omega, obj.i_vel, obj.b_tau}, ...
                {ddX}, ...
                {'Rotation_matrix','Omega_body','Vel_inertial','Wrench_body'}, ...
                {'State_derivative'});
        end
    end

end