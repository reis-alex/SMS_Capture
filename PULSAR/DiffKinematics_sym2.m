function [Bij,Bi0,P0,pm] = DiffKinematics_sym2(R0, r0, rL, e, g, robot)

import casadi.*
SkewSym = @(x)[0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0];

% Number of links
n = robot.n_links_joints;

% Allocate
for i = 1:n
    for j = 1:n
        Bij{i}{j} = SX.zeros(6,6);
    end
end
Bi0 = SX.sym('Bi0',6,6,n);

% --------------------------------------------
% Convert link positions to BASE frame
% --------------------------------------------
rL_base = SX.zeros(3,n);
for i = 1:n
    rL_base(:,i) = R0' * ( rL(1:3,i) - r0 );  % relative position, expressed in base frame
end

% --------------------------------------------
% Build Bij in base frame
% Bij relates link i to link j
% --------------------------------------------
for j = 1:n
    for i = 1:n
        if robot.con.branch(i,j) == 1
            Bij{i}{j}(1:6,1:6) = ...
                [ SX.eye(3), SX.zeros(3,3);
                  SkewSym( rL_base(:,j) - rL_base(:,i) ), SX.eye(3) ];
        else
            Bij{i}{j}(1:6,1:6) = SX.zeros(6,6);
        end
    end
end

% --------------------------------------------
% Build Bi0 (link i relative to base)
% --------------------------------------------
for i = 1:n
    Bi0{i}(:,:) = ...
        [ SX.eye(3), SX.zeros(3,3);
          SkewSym( rL_base(:,i) ), SX.eye(3) ];
end

% --------------------------------------------
% Base twist propagation matrix
% --------------------------------------------
P0 = [R0, SX.zeros(3,3);
      SX.zeros(3,3), SX.eye(3)];

% --------------------------------------------
% pm vectors
% --------------------------------------------
pm = SX.zeros(6,n);
for i = 1:n
    if robot.joints(i).type == 1     % revolute
        pm(:,i) = [e(:,i); cross(e(:,i), g(:,i))];
    elseif robot.joints(i).type == 2 % prismatic
        pm(:,i) = [SX.zeros(3,1); e(:,i)];
    else                              % fixed
        pm(:,i) = SX.zeros(6,1);
    end
end

end
