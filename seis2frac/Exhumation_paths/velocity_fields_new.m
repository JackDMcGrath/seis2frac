function [vx,vy,vz,AF]=velocity_fields_new(bearing,origin,xx,yy,zz,strikeSlip,converge,dipSlip,exhum,af)
% New method - constant convergence.
% Everything on the other side of the AF set to 0

dip=atand((af(1,3)-af(3,3))/(af(1,2)-af(3,2)));

vx=ones(size(xx))*strikeSlip;
vy=ones(size(yy))*converge;

[exrot]=point_rotation(exhum(:,1:2),90-bearing,'strike');
exhum(:,1:2)=exrot(:,1:2);
exhum_interp=scatteredInterpolant(exhum(:,1),exhum(:,2),exhum(:,3));
exhum_interp.ExtrapolationMethod='nearest'; % options nearest or linear
ex_vz=exhum_interp(xx(:),yy(:))*1e-6;
vz=reshape(ex_vz,size(xx));

%% Identify the Australian plate
% Find the normal of the AF
%https://au.mathworks.com/matlabcentral/answers/75952-position-of-points-relative-to-a-plane

N=null(-af);
planeOrigin=[0;0;0];
points=[xx(:),yy(:),zz(:)];

aus=reshape((dot(points'-planeOrigin, repmat(N,[1,size(points,1)]))>0),size(yy));
vx(aus)=0;
vy(aus)=0;
vz(aus)=0;
vx=reshape(vx,size(yy));
vy=reshape(vy,size(yy));
vz=reshape(vz,size(yy));

%% Identify AF region, and assign that a velocity

faultWidth=1; % Thickness of the AF region

afvel_min=(aus==0);

planeOrigin=[0;-(faultWidth/sind(dip));0];

afvel_max=reshape((dot(points'-planeOrigin, repmat(N,[1,size(points,1)]))<0),size(yy));

%overlap of afvels are the points in the AF region
afvel=(afvel_min~=afvel_max);

vy(afvel)=cosd(dip)*dipSlip;
vz(afvel)=sind(dip)*dipSlip;
vy=reshape(vy,size(yy));
vz=reshape(vz,size(yy));

% Coordinates of the boundary of the AF region
AF.surface=af;
af(:,2)=af(:,2)+planeOrigin(2);

% Structure with AF region details
AF.surface=[AF.surface;af];
AF.normal=N;
AF.vy=cosd(dip)*dipSlip;
AF.vz=sind(dip)*dipSlip;
AF.dip=dip;
end