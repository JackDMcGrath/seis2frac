function [thickness]=depth_surf(orig,seis_depths,bearing,dip,plot_fig);
%function plots upper and lower seismicity as a surface, then works out the
%thickness of the seismic zone for a rock volume travelling not-vertically
%through the seismogenic zone.

% plot_fig=0;

x=seis_depths(~isnan(seis_depths(:,1)),1);
y=seis_depths(~isnan(seis_depths(:,1)),2);
faces=delaunay(x,y);

z_up=seis_depths(~isnan(seis_depths(:,1)),3);
z_down=seis_depths(~isnan(seis_depths(:,1)),4);

v_up=[x,y,z_up];
v_down=[x,y,z_down];

vu1=v_up(faces(:,1),:);
vu2=v_up(faces(:,2),:);
vu3=v_up(faces(:,3),:);

vd1=v_down(faces(:,1),:);
vd2=v_down(faces(:,2),:);
vd3=v_down(faces(:,3),:);

dir=[-cosd(-bearing)*(50/tand(dip)) -sind(-bearing)*(50/tand(dip)) 50]; % Direction from bin center, gaining 50km height

intersectu = TriangleRayIntersection(orig, dir, vu1, vu2, vu3);
intersectd = TriangleRayIntersection(orig, dir, vd1, vd2, vd3);
% if isempty(intersectu)
% intersectu = TriangleRayIntersection(orig, -dir, vu1, vu2, vu3);
% end
% if isempty(intersectd)
%     intersectd = TriangleRayIntersection(orig, -dir, vd1, vd2, vd3);
% end

if isempty(find(intersectu)) || isempty(find(intersectd))
    intersectu = TriangleRayIntersection(orig-dir, dir, vu1, vu2, vu3);
    intersectd = TriangleRayIntersection(orig-dir, dir, vd1, vd2, vd3);
end


if plot_fig==1
    figure(); clf;
    trisurf(faces,x,y,z_up, intersectu*1.0,'FaceAlpha', 0.9)
    hold on;
    trisurf(faces,x,y,z_down, intersectd*1.0,'FaceAlpha', 0.9)
    line('XData',orig(1)+[0 dir(1)],'YData',orig(2)+[0 dir(2)],'ZData',...
        orig(3)+[0 dir(3)],'Color','r','LineWidth',3)
    set(gca, 'CameraPosition', [106.2478  -35.9079  136.4875])
    %set(gco,'EdgeColor','none');
end

upper=mean([vu1(find(intersectu),3);vu2(find(intersectu),3);vu3(find(intersectu),3)]);
lower=mean([vd1(find(intersectd),3);vd2(find(intersectd),3);vu3(find(intersectd),3)]);
thickness=upper-lower;


if isnan(thickness) % Try replacing any that have just missed the surface with a vertical thickness
    thickness=griddata(seis_depths(~isnan(seis_depths(:,1)),1),seis_depths(~isnan(seis_depths(:,1)),2),seis_depths(~isnan(seis_depths(:,1)),3),orig(1),orig(2))- ...
        griddata(seis_depths(~isnan(seis_depths(:,1)),1),seis_depths(~isnan(seis_depths(:,1)),2),seis_depths(~isnan(seis_depths(:,1)),4),orig(1),orig(2));
%     if isnan(thickness)
%         line('XData',orig(1)+[0 dir(1)],'YData',orig(2)+[0 dir(2)],'ZData',...
%             orig(3)+[0 dir(3)],'Color','r','LineWidth',3)
%         line('XData',orig(1)-[0 dir(1)],'YData',orig(2)-[0 dir(2)],'ZData',...
%             orig(3)-[0 dir(3)],'Color','r','LineWidth',3)
%     else
%         line('XData',orig(1)+[0 dir(1)],'YData',orig(2)+[0 dir(2)],'ZData',...
%             orig(3)+[0 dir(3)],'Color','b','LineWidth',3)
%         line('XData',orig(1)-[0 dir(1)],'YData',orig(2)-[0 dir(2)],'ZData',...
%             orig(3)-[0 dir(3)],'Color','b','LineWidth',3)
%     end
%     
% else
%     line('XData',orig(1)+[0 dir(1)],'YData',orig(2)+[0 dir(2)],'ZData',...
%         orig(3)+[0 dir(3)],'Color','g','LineWidth',3)
%     line('XData',orig(1)-[0 dir(1)],'YData',orig(2)-[0 dir(2)],'ZData',...
%         orig(3)-[0 dir(3)],'Color','g','LineWidth',3)
end
plot3(orig(1),orig(2),orig(3),'k.','MarkerSize',10)
end