function [display]=update_display_robot(display,robot,R0, r0, qm)
%--- Kinematics ---%
%Kinematics
[RJ,RL,rJ,rL,e,g]=Kinematics(R0,r0,qm,robot);
T0=[R0,r0;zeros(1,3),1];

for ind=1:robot.n_links_joints
    if ~isempty(robot.links(ind).visual.vertices)
        ind_parent=(robot.links(ind).parent_joint);
        T_joint=[RL(:,:,ind_parent) rJ(:,ind_parent); 0 0 0 1];
        T=T_joint*robot.links(ind).visual.T;
        v=[robot.links(ind).visual.vertices ones(length(robot.links(ind).visual.vertices),1)];
        vert=T*v';
        
        display.body(ind).Vertices=vert(1:3,:)';
    end
%     display.COM(ind).XData=rL(1,ind);
%     display.COM(ind).YData=rL(2,ind);
%     display.COM(ind).ZData=rL(3,ind);
%     display.Joint(ind).XData=rJ(1,ind);
%     display.Joint(ind).YData=rJ(2,ind);
%     display.Joint(ind).ZData=rJ(3,ind);
end
ind=ind+1;
if ~isempty(robot.base_link.visual.vertices)
    v=[robot.base_link.visual.vertices ones(length(robot.base_link.visual.vertices),1)];
    vert=T0*robot.base_link.visual.T*v';
    display.body(ind).Vertices=vert(1:3,:)';
end
drawnow();