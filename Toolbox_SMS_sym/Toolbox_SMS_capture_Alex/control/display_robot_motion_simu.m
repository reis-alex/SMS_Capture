function [display]=display_robot_motion_simu(robot,log,step,t_start)
% Display a SPART robot model created from a URDF file.
w=get(0,'ScreenSize');

figure('Position',[0 0 w(3)*2/3 w(4)*2/3]),
% figure()
hold on;

%--- Kinematics ---%
%Kinematics

R0 = log.getElement('R0');
r0 = log.getElement('r0');
qm = log.getElement('q_simu');
rL_target = log.getElement('rL_target');
rL_EE = log.getElement('rL_EE');
u_EE = log.getElement('u_EE');
u_target = log.getElement('u_t');

% [~,RL,rJ,rL,~,~]=Kinematics(R0.Values.Data(:,:,1),r0.Values.Data(1,:)',qm.Values.Data(1,:),robot);

[~,RL,rJ,rL,~,~]=Kinematics(R0.Values.Data(:,:,t_start),r0.Values.Data(t_start,:)',qm.Values.Data(t_start,:),robot);


for ind=1:robot.n_links_joints
    if ~isempty(robot.links(ind).visual.vertices)
        ind_parent=(robot.links(ind).parent_joint);
        T_joint=[RL(:,:,ind_parent) rJ(:,ind_parent); 0 0 0 1];
        T=T_joint*robot.links(ind).visual.T;
        v=[robot.links(ind).visual.vertices ones(length(robot.links(ind).visual.vertices),1)];
        vert=T*v';
        if strcmp(robot.links(ind).stif.type,'flexible')
            FaceAlpha=0.4;
        else
            FaceAlpha=0.9;
        end
        
        display.body(ind)=patch('Faces',robot.links(ind).visual.faces,...
            'Vertices',vert(1:3,:)',...
            'FaceColor',       [0.8 0.8 1.0], ...
            'EdgeColor',       'none',        ...
            'FaceLighting',    'gouraud',     ...
            'AmbientStrength', 0.15, ...
            'FaceAlpha',FaceAlpha);
    end
    display.COM(ind)=plot3(rL(1,ind),rL(2,ind),rL(3,ind),'LineWidth',2,'Color','k','Marker','+');
    display.Joint(ind)=plot3(rJ(1,ind),rJ(2,ind),rJ(3,ind),'LineWidth',2,'Color','r','Marker','o');
    if ind==1
        camlight('headlight');
    end
    material('dull');
    axis('image');
    view([-135 35]);
end
ind=ind+1;
if ~isempty(robot.base_link.visual.vertices)
    v=[robot.base_link.visual.vertices ones(length(robot.base_link.visual.vertices),1)];
    vert=robot.base_link.visual.T*v';
    display.body(ind)=patch('Faces',robot.base_link.visual.faces,...
        'Vertices',vert(1:3,:)',...
        'FaceColor',       [0.9 0.4 0.4], ...
        'EdgeColor',       'none',        ...
        'FaceLighting',    'gouraud',     ...
        'AmbientStrength', 0.15, ...
        'FaceAlpha',.8);
end




% for t=1:step:length(R0.Values.Data)
for t=t_start:step:length(R0.Values.Data)
    hold on;
    [~,RL,rJ,rL,~,~]=Kinematics(R0.Values.Data(:,:,t),r0.Values.Data(t,:)',qm.Values.Data(t,:),robot);
    
    color = [1 1 t/length(R0.Values.Data)];

    for ind=4:robot.n_links_joints
        if ~isempty(robot.links(ind).visual.vertices)
            ind_parent=(robot.links(ind).parent_joint);
            T_joint=[RL(:,:,ind_parent) rJ(:,ind_parent); 0 0 0 1];
            T=T_joint*robot.links(ind).visual.T;
            v=[robot.links(ind).visual.vertices ones(length(robot.links(ind).visual.vertices),1)];
            vert=T*v';
            if strcmp(robot.links(ind).stif.type,'flexible')
                FaceAlpha=0.4;
            else
                FaceAlpha=0.9;
            end
            
            display.body(ind)=patch('Faces',robot.links(ind).visual.faces,...
                'Vertices',vert(1:3,:)',...
                'FaceColor',       color, ...
                'EdgeColor',       'none',        ...
                'FaceLighting',    'gouraud',     ...
                'AmbientStrength', 0.15, ...
                'FaceAlpha',FaceAlpha);
        end
        display.COM(ind)=plot3(rL(1,ind),rL(2,ind),rL(3,ind),'LineWidth',2,'Color',color,'Marker','+');
        display.Joint(ind)=plot3(rJ(1,ind),rJ(2,ind),rJ(3,ind),'LineWidth',2,'Color',color,'Marker','o');
        if ind==1
            camlight('headlight');
        end
        material('dull');
        axis('image');
        view([-135 35]);
    end
    ind=ind+1;
    if ~isempty(robot.base_link.visual.vertices)
        % v=[robot.base_link.visual.vertices ones(length(robot.base_link.visual.vertices),1)];
        % vert=robot.base_link.visual.T*v';

        T_0=[R0.Values.Data(:,:,t) r0.Values.Data(t,:)'; 0 0 0 1];
        T=T_0*robot.base_link.visual.T;
        v=[robot.base_link.visual.vertices ones(length(robot.base_link.visual.vertices),1)];
        vert=T*v';

        display.body(ind)=patch('Faces',robot.base_link.visual.faces,...
            'Vertices',vert(1:3,:)',...
            'FaceColor',       color, ...
            'EdgeColor',       'none',        ...
            'FaceLighting',    'gouraud',     ...
            'AmbientStrength', 0.55, ...
            'FaceAlpha',.8);
    end

end

hold on;
% colormap(jet);
% % scatter3(xee(1,:),xee(2,:),xee(3,:),4,indice,'filled');
% scatter3(rL_EE.Values.Data(:,1),rL_EE.Values.Data(:,2),rL_EE.Values.Data(:,3),15, sqrt(sum(u_EE.Values.Data.^2,2)),'filled');
% % colorbar;
% 
% hold on;
% colormap(jet);
% scatter3(rL_target.Values.Data(:,1),rL_target.Values.Data(:,2),rL_target.Values.Data(:,3),15, sqrt(sum(u_target.Values.Data.^2,2)),'filled');
% colorbar;

f = gcf;

ax1 = f.CurrentAxes;
% axis equal
% ax2 = axes;


% hLink = linkprop([ax1,ax2],{'X','Y','Z','X','Y','Z'});

% ax2.Visible = 'off';
% ax2.XTick = [];
% ax2.YTick = [];

colormap(ax1,'jet');
% colormap(ax2,'summer');

scatter3(ax1,rL_EE.Values.Data(:,1),rL_EE.Values.Data(:,2),rL_EE.Values.Data(:,3),15, sqrt(sum(u_EE.Values.Data.^2,2)),'filled');
hold on;
% scatter3(ax2,rL_target.Values.Data(:,1),rL_target.Values.Data(:,2),rL_target.Values.Data(:,3),15, sqrt(sum(u_target.Values.Data.^2,2)),'filled');
cb1 = colorbar(ax1,'Position',[0.1 0.1 0.05 0.815]);
% cb2 = colorbar(ax2,'Position',[0.9 0.1 0.05 0.815]);

cb1.Label.String = 'norm uEE';
% cb2.Label.String = 'norm ut';
cb1.Label.FontSize = 14;
% cb2.Label.FontSize = 14;


xlabel('X(m)')
ylabel('Y(m)')
zlabel('Z(m)')
%title(tl);
grid on