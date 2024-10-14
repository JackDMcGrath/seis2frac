function [vx,vy,vz,AF]=updateVelocityField(xx,yy,zz,vx,vy,vz,AF,width,height,compress,xyz)
% Function to update the velocity fields so that when a rock package hits
% the Alpine Fault, subsequent rocks are forced up at the same rate earlier
% planOrigin=impact point
% width=horizontal seperation of rays
% height=vertical seperation of rays
% Compression= width/compression is how much you will allow rocks to
% compress as they hit the fault

points=[xx(:),yy(:),zz(:)];
N=AF.normal;

newHeight=height/compress; % This is how far out we push the AF region
planeOrigin=xyz;
planeOrigin(2)=planeOrigin(2)-(newHeight/sind(AF.dip)); % New plane would pass through here


aus=reshape((dot(points'-[0;0;0], repmat(N,[1,size(points,1)]))>0),size(yy));
afvel_max=reshape((dot(points'-planeOrigin, repmat(N,[1,size(points,1)]))<0),size(yy));
afvel_left=reshape(points(:,1)<(planeOrigin(1)-0.5*width),size(yy));
afvel_right=reshape(points(:,1)>(planeOrigin(1)+0.5*width),size(yy));
afvel_base=reshape(points(:,3)<(planeOrigin(3)),size(yy));

% Where these all overlap (ie =0), this is the new area 
new=(afvel_right+afvel_left+aus+afvel_max+afvel_base)==0;

vy(new)=AF.vy;
vz(new)=AF.vz;

% Find coordinates of new region volume [TL,TR,BR,BL]
% New protruding section
TL=[planeOrigin(1)-0.5*width, planeOrigin(2)-planeOrigin(3)/tand(AF.dip), 0];
TR=[planeOrigin(1)+0.5*width, planeOrigin(2)-planeOrigin(3)/tand(AF.dip), 0];
BR=[planeOrigin(1)+0.5*width, planeOrigin(2), planeOrigin(3)];
BL=[planeOrigin(1)-0.5*width, planeOrigin(2), planeOrigin(3)];
surface=[TL;TR;BR;BL];

% Section joining to the pre-existing region (ie on fault (the weld...))
surface2=surface;
surface2(:,2)=surface2(:,2)+(newHeight/sind(AF.dip));

new_volume=[surface;surface2];

faces=boundary(new_volume(:,1),new_volume(:,2),new_volume(:,3));
faces=faces+size(AF.surface,1);

AF.surface=[AF.surface;new_volume];
AF.faces=[AF.faces;faces];
end

