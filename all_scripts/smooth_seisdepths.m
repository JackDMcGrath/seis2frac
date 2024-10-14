function [seis_depths]=smooth_seisdepths(seis_depths);

% Function to take the quadtree seisdepths, and to interpolate to more
% regular grid

x=seis_depths(~isnan(seis_depths(:,1)),1);
y=seis_depths(~isnan(seis_depths(:,1)),2);
faces=delaunay(x,y);

for ii=1:length(faces)
    dist(ii)=sqrt((seis_depths(faces(ii,1),1)-seis_depths(faces(ii,2),1))^2+(seis_depths(faces(ii,1),2)-seis_depths(faces(ii,2),2))^2);
        dist(ii+1)=sqrt((seis_depths(faces(ii,2),1)-seis_depths(faces(ii,3),1))^2+(seis_depths(faces(ii,2),2)-seis_depths(faces(ii,3),2))^2);
            dist(ii+2)=sqrt((seis_depths(faces(ii,1),1)-seis_depths(faces(ii,2),3))^2+(seis_depths(faces(ii,1),2)-seis_depths(faces(ii,3),2))^2);
end

grid_spacing=floor(mean(dist)/5)*5; % set grid to mean spacing of data nodes (to lowest 5)

ylim=[-80:grid_spacing:10];
xlim=[0:grid_spacing:280];

[xx,yy]=meshgrid(xlim,ylim);

upper=scatteredInterpolant(x,y,seis_depths(:,3));
upper.Method='natural';
zup=upper(xx(:),yy(:));
lower=scatteredInterpolant(x,y,seis_depths(:,4));
zlow=lower(xx(:),yy(:));
seis_depths=[xx(:),yy(:),zup,zlow];
end