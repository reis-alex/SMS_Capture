
function [mtot] =mass_tot(robot)

mass_total=robot.base_link.mass;

%Add contributions of links
for i=1:robot.n_links_joints
    mass_total=mass_total+robot.links(i).mass;
end

mtot =mass_total;