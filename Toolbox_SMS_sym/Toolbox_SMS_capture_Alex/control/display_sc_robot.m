function [display_rb,display_sc]=display_sc_robot(robot,sc,R_rs2rt, rt, qm)
% Display a SPART robot model created from a URDF file.
display_rb=display_robot(robot,R_rs2rt,rt,qm);
%--- Kinematics ---%
%Kinematics
[RJ,RL,rJ,rL,e,g]=Kinematics(eye(3),[0 0 0]',zeros(sc.n_q,1),sc);


for ind=1:sc.n_links_joints
    if ~isempty(sc.links(ind).visual.vertices)
        ind_parent=(sc.links(ind).parent_joint);
        T_joint=[RL(:,:,ind_parent) rJ(:,ind_parent); 0 0 0 1];
        T=T_joint*sc.links(ind).visual.T;
        v=[sc.links(ind).visual.vertices ones(length(sc.links(ind).visual.vertices),1)];
        vert=T*v';
        if strcmp(sc.links(ind).stif.type,'flexible')
            FaceAlpha=0.4;
        else
            FaceAlpha=0.9;
        end
        
        display_sc.body(ind)=patch('Faces',sc.links(ind).visual.faces,...
            'Vertices',vert(1:3,:)',...
            'FaceColor',       [0.8 0.8 1.0], ...
            'EdgeColor',       'none',        ...
            'FaceLighting',    'gouraud',     ...
            'AmbientStrength', 0.15, ...
            'FaceAlpha',FaceAlpha);
    end
    display_sc.COM(ind)=plot3(rL(1,ind),rL(2,ind),rL(3,ind),'LineWidth',2,'Color','k','Marker','+');
    display_sc.Joint(ind)=plot3(rJ(1,ind),rJ(2,ind),rJ(3,ind),'LineWidth',2,'Color','r','Marker','o');
    if ind==1
        camlight('headlight');
    end
    material('dull');
    axis('image');
    view([-135 35]);
end
ind=ind+1;
if ~isempty(sc.base_link.visual.vertices)
    v=[sc.base_link.visual.vertices ones(length(sc.base_link.visual.vertices),1)];
    vert=sc.base_link.visual.T*v';
    display_sc.body(ind)=patch('Faces',sc.base_link.visual.faces,...
        'Vertices',vert(1:3,:)',...
        'FaceColor',       [0.9 0.4 0.4], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15, ...
        'FaceAlpha',.8);
end

xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
grid on