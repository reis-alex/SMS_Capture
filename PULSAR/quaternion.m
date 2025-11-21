classdef quaternion
    properties

    end
    methods
        function [qplus, R] = integrate(obj, omega_body, q0, dt)
            % Integrates quaternion from body angular velocity.
            %
            % Inputs:
            %   q0          : initial quaternion [q0; q1; q2; q3] (unit norm)
            %   dt          : integration step
            %
            % Outputs:
            %   q : updated quaternion
            %   R : rotation matrix

            wx = omega_body(1);
            wy = omega_body(2);
            wz = omega_body(3);
            
            Omega = [  0,   -wx, -wy, -wz;
                        wx,   0,    wz, -wy;
                        wy,  -wz,   0,   wx;
                        wz,   wy,  -wx,   0 ];

            qplus = q0 + 0.5 * Omega * q0 * dt;
            qplus = qplus / norm(qplus);  % renormalize

            % Compute rotation matrices
            R = zeros(3,3);

            q0 = qplus(1);
            q1 = qplus(2); 
            q2 = qplus(3); 
            q3 = qplus(4);

            R = [1-2*(q2^2+q3^2),   2*(q1*q2 - q0*q3),   2*(q1*q3 + q0*q2);
                2*(q1*q2 + q0*q3), 1-2*(q1^2+q3^2),     2*(q2*q3 - q0*q1);
                2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1),   1-2*(q1^2+q2^2)];
        end

        function q = euler_to_quat(obj, phi, theta, psi)
            % Convert Euler angles (rad) → quaternion [q0; q1; q2; q3]
            % Convention: 3-2-1 (yaw-pitch-roll)
            c1 = cos(psi/2);  s1 = sin(psi/2);
            c2 = cos(theta/2); s2 = sin(theta/2);
            c3 = cos(phi/2);  s3 = sin(phi/2);

            q = [c1*c2*c3 + s1*s2*s3;
                s3*c2*c1 - c3*s2*s1;
                c3*s2*c1 + s3*c2*s1;
                c3*c2*s1 - s3*s2*c1 ];

            q = q / norm(q); % ensure normalization
        end


        function q = rotm_to_quat(obj,R)
            % Convert rotation matrix → quaternion [q0; q1; q2; q3]
            tr = trace(R);
            if tr > 0
                S = sqrt(tr + 1.0) * 2;
                q0 = 0.25 * S;
                q1 = (R(3,2) - R(2,3)) / S;
                q2 = (R(1,3) - R(3,1)) / S;
                q3 = (R(2,1) - R(1,2)) / S;
            else
                % Handle cases when trace is small (to avoid division by zero)
                [~, i] = max([R(1,1), R(2,2), R(3,3)]);
                switch i
                    case 1
                        S = sqrt(1.0 + R(1,1) - R(2,2) - R(3,3)) * 2;
                        q0 = (R(3,2) - R(2,3)) / S;
                        q1 = 0.25 * S;
                        q2 = (R(1,2) + R(2,1)) / S;
                        q3 = (R(1,3) + R(3,1)) / S;
                    case 2
                        S = sqrt(1.0 + R(2,2) - R(1,1) - R(3,3)) * 2;
                        q0 = (R(1,3) - R(3,1)) / S;
                        q1 = (R(1,2) + R(2,1)) / S;
                        q2 = 0.25 * S;
                        q3 = (R(2,3) + R(3,2)) / S;
                    case 3
                        S = sqrt(1.0 + R(3,3) - R(1,1) - R(2,2)) * 2;
                        q0 = (R(2,1) - R(1,2)) / S;
                        q1 = (R(1,3) + R(3,1)) / S;
                        q2 = (R(2,3) + R(3,2)) / S;
                        q3 = 0.25 * S;
                end
            end
            q = [q0; q1; q2; q3];
            q = q / norm(q);
        end


        function angles = quat_to_angles(obj, q)
            % q = [q0; q1; q2; q3] (scalar first)
            % 3-2-1 Euler angles: phi-roll, theta-pitch, psi-yaw

            q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

            angles(1) = atan2(2*(q0*q1 + q2*q3), 1 - 2*(q1^2 + q2^2));  % roll
            angles(2) = asin(2*(q0*q2 - q3*q1));                        % pitch
            angles(3) = atan2(2*(q0*q3 + q1*q2), 1 - 2*(q2^2 + q3^2));  % yaw
        end
    end
end