function [p1,p2,p3,p4]=find_corners2(center,dim,vector)
% Build square around centre point, and define corners
% p1,p2,p3,p4 are TL, TR, BR, BL respectively

xyvec=vector(:,[1 2]);
xzvec=vector(:,[1 3]);
yzvec=vector(:,[2 3]);

%% X-Y Plane

V=sqrt(xyvec(:,1).^2+xyvec(:,2).^2);
dip=atand(xyvec(:,2)./xyvec(:,1));
dip23=dip-atand(dim/V);
dip14=dip+atand(dim/V);
V2=sqrt(V.^2+dim^2);
v23x=cosd(dip23).*V2;
v23y=sind(dip23).*V2;
v14x=cosd(dip14).*V2;
v14y=sind(dip14).*V2;

%% X-Z Plane

V=sqrt(xzvec(:,1).^2+xzvec(:,2).^2);
dip=atand(xzvec(:,2)./xzvec(:,1));
dip23=dip-atand(dim/V);
dip14=dip+atand(dim/V);
V2=sqrt(V.^2+dim^2);
v23x=cosd(dip23).*V2;
v23z=sind(dip23).*V2;
v14x=cosd(dip14).*V2;
v14z=sind(dip14).*V2;

%% Y-Z Plane

V=sqrt(yzvec(:,1).^2+yzvec(:,2).^2);
dip=atand(yzvec(:,2)./yzvec(:,1));
dip23=dip-atand(dim/V);
dip14=dip+atand(dim/V);
V2=sqrt(V.^2+dim^2);
v23y=cosd(dip23).*V2;
v23z=sind(dip23).*V2;
v14y=cosd(dip14).*V2;
v14z=sind(dip14).*V2;



end