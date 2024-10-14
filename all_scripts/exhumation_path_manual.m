dbstop if error
warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
%% Prepare Velocity Fields
bearing=56.8962;
origin=[168.685,-44.115,0];
fault=[0,0;324.2663,0];

lon=-50:10:250;
lat=-100:10:20;
[xx,yy]=meshgrid(lon,lat);

gps=readmatrix('gps.csv');
gps=gps(1:937,:);
[in]=inpolygon(gps(:,1),gps(:,2),[169,169,173.5,173.5,169],[-45,-42,-42,-45,-45]);
gps=gps(in,:);
gps(:,3:4)=[gps(:,3)*cosd(90-bearing)+gps(:,4)*sind(90-bearing),-gps(:,3)*sind(90-bearing)-gps(:,4)*cosd(90-bearing)];
gps(:,1:2)=ll2utm(gps(:,(1:2)),origin(:,1:2));
[gpsrot]=point_rotation(gps(:,1:2),90-bearing,'strike');
gps(:,1:2)=gpsrot(:,1:2);
gpsref(1)=griddata(gps(:,1),gps(:,2),gps(:,3),100,0);
gpsref(2)=griddata(gps(:,1),gps(:,2),gps(:,4),100,0);
% gps(:,3:4)=gps(:,3:4)-gpsref;

exhum=readmatrix('C:\Jacko\NZ 2020\seis2frac3\Modeled_exhumation_rates_grid_points.csv');
exhum([11,1:10],3)=exhum(12:22,3); %Make AUS nodes equal PACmax
exhum(:,(1:2))=ll2utm(exhum(:,(1:2)),origin(:,1:2));
[exrot]=point_rotation(exhum(:,1:2),90-bearing,'strike');
exhum(:,1:2)=exrot(:,1:2);

dip_slip=readmatrix('C:\Jacko\NZ 2020\seis2frac3\dip_slip_rates.csv');
dip_slip(:,(2:3))=ll2utm(dip_slip(:,(2:3)),origin(:,1:2));
[diprot]=point_rotation(dip_slip(:,2:3),90-bearing,'strike');
dip_slip(:,2:3)=diprot(:,1:2);
dip_slip(:,3)=0; % Put all slip rates on the "fault"
sliprate=dip_slip(~isnan(dip_slip(:,4)),[2 3 4]); % Slip rate along the fault
dipv=dip_slip(~isnan(dip_slip(:,5)),[2 3 5]);
dipn=dipv;
dipv(:,3)=sind(60)*dipv(:,3); % Uplift due to fault slip
dipn(:,3)=cos(60)*dipn(:,3); % Convergence due to fault slip
slip_distance=60; % Distance at which fault velocity contribution is removed

sliprate_grid=zeros(size(xx));
dipv_grid=zeros(size(xx));
dipn_grid=zeros(size(xx));


f_ix=find(yy(:,1)==0); % Find row equating to fault
sliprate_grid(f_ix,:)=interp1(sliprate(:,1),sliprate(:,3),xx(f_ix,:)); % Linear interp to grid
[mdl]=fitlm(xx(1,:),sliprate_grid(f_ix,:)); % Linear model though these (given overall trend is vaguely linear
sliprate_grid(f_ix,:)=mdl.Coefficients.Estimate(1)+xx(1,:)*mdl.Coefficients.Estimate(2);
dipv_grid(f_ix,:)=interp1(dipv(:,1),dipv(:,3),xx(f_ix,:),'nearest');
dipn_grid(f_ix,:)=interp1(dipn(:,1),dipn(:,3),xx(f_ix,:),'nearest');

for ii=1:size(xx,2)
    sliprate_grid(:,ii)=(2*sliprate_grid(f_ix,ii)/pi)*atan(abs(yy(:,ii))*1e3/15e3);
    dipv_grid(:,ii)=(2*dipv_grid(f_ix,ii)/pi)*atan(yy(:,ii)*1e3/15e3);
    dipn_grid(:,ii)=(2*dipn_grid(f_ix,ii)/pi)*atan(yy(:,ii)*1e3/15e3);
end

sliprate_grid=sliprate_grid-mdl.Coefficients.Estimate(1)+xx(1,:)*mdl.Coefficients.Estimate(2);
dipv_grid=dipv_grid+interp1(dipv(:,1),dipv(:,3),xx(f_ix,:),'nearest');
dipn_grid=dipn_grid+interp1(dipn(:,1),dipn(:,3),xx(f_ix,:),'nearest');


clear gpsrot exrot diprot


%% Prepare depth surfaces

% load depth_file.mat
% % [seisrot]=point_rotation(seis_depths(:,1:2),90-bearing,'strike'); % DATA ALREADY ROTATED
% % seis_depths(:,1:2)=seisrot(:,1:2);
% % seis_depths(isnan(seis_depths(:,3)),:)=[];
% ix=find(~isnan(seis_depths(:,3)));
% seis_depths(:,3)=griddata(seis_depths(ix,1),seis_depths(ix,2),seis_depths(ix,3),seis_depths(:,1),seis_depths(:,2),'cubic');
% 
% [seis_depths]=smooth_seisdepths(seis_depths);

load depth_surfaceTmean.mat;

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

%% pt = [ Lon Lat Depth uppersurface_depth uppertime lowersurface_depth lowersurface_time]
pt=[142.0609,-7.4839,0,0,0,0,0]; % Crawford's Knob
% pt=[141.6346,-2.6302,0,0,0,0,0]; % Franz Josef Planar
dir=[0,0,0]; % Vector - initially surface normal
pt=[140, 0,0,0,0,0,0]; % Crawford's Knob

poly_width=5; %width of poly in km
p=nan(4,3);
[p(1,:),p(2,:),p(3,:),p(4,:)]=find_corners(pt(1,1:3),poly_width);

depth=pt(end,3);
max_depth=30;

timestep=1e5;
step=1;

f= figure;
hold on
xlabel('X');ylabel('Y');zlabel('Z')
plot3(pt(end,1),pt(end,2),pt(end,3),'r.');
view([-45,30])
grid

searchupper=1;
searchlower=0;
end_search=0;

while depth<max_depth
    
    [dx]=griddata(gps(:,1),gps(:,2),gps(:,3),pt(end,1),pt(end,2));
    [dy]=griddata(gps(:,1),gps(:,2),gps(:,4),pt(end,1),pt(end,2));
    [dz]=griddata(exhum(:,1),exhum(:,2),exhum(:,3),pt(end,1),pt(end,2));
    [dx_slip]=griddata(xx(:),yy(:),sliprate_grid(:),pt(end,1),pt(end,2));
    [dy_slip]=griddata(xx(:),yy(:),dipn_grid(:),pt(end,1),pt(end,2));
    [dz_slip]=griddata(xx(:),yy(:),dipv_grid(:),pt(end,1),pt(end,2));
    if isnan(dz)
        dz=0.1;
    end
    dir(step+1,:)=-[dx*1e-6*timestep,dy*1e-6*timestep,dz*1e-6*timestep];
%     dir(step+1,:)=-[dx*1e-6*timestep,dy*1e-6*timestep,(dz+dz_slip)*1e-6*timestep]; % Using the dipslip
%     dir(step+1,:)=-[(dx+dx_slip)*1e-6*timestep,(dy+dy_slip)*1e-6*timestep,dz*1e-6*timestep]; % Using the sliprates
%     dir(step+1,:)=-[(dx+dx_slip)*1e-6*timestep,(dy+dy_slip)*1e-6*timestep,(dz+dz_slip)*1e-6*timestep]; % Using the sliprates + dipslip
    intersectusurf=zeros(size(vu1,1),1);
    if searchupper==1
        intersectu = TriangleRayIntersection(pt(end,1:3), dir(step+1,:), vu1, vu2, vu3,'LineType','Segment'); % Linetype makes backcheck redundant
        if ~isempty(find(intersectu))
%             backcheck = TriangleRayIntersection(pt(end,1:3)+dir(step+1,:), -dir(step+1,:), vu1, vu2, vu3); % backcheck needed to ensure we've passed the surface
%             if ~isempty(find(backcheck)) % If the return vector also passes the surface
%                 triangles=find(intersectu); % Find intersection of forward vector
%                 back_triangles=find(backcheck); % Find intersetion of back vector
%                 intersection=find(ismember(triangles,back_triangles));  % Search for the same intersection in both directions
%                 if ~isempty(intersection) % If you have passed through the layer, this will not be empty
%                     ix=triangles(intersection); % Recover ID of intersection triangle
%                     ix=faces(ix,:); % Get corners of said triangle
%                     pt(end,4)=mean([seis_depths(ix(1),3),seis_depths(ix(2),3),seis_depths(ix(3),3)]); % change to weighted mean
%                     pt(end,5)=step;
                    searchupper=0;
                    searchlower=1;
%                     intersectu=zeros(size(intersectu));
%                     intersectu(triangles(ismember(triangles,find(backcheck))))=1;
intersectusurf(find(intersectu))=1;
                    trisurf(faces,x,y,z_up, intersectusurf*1.0,'FaceAlpha', 0.9);
%                 end
%             end
        end
    end
    
    if searchlower==1
        intersectd = TriangleRayIntersection(pt(end,1:3), dir(step+1,:), vd1, vd2, vd3);
        if ~isempty(find(intersectd))
            backcheck = TriangleRayIntersection(pt(end,1:3)+dir(step+1,:), -dir(step+1,:), vd1, vd2, vd3);
            if ~isempty(find(backcheck)) % If the return vector also passes the surface
                triangles=find(intersectd); % Find intersection of forward vector
                back_triangles=find(backcheck); % Find intersetion of back vector
                intersection=find(ismember(triangles,back_triangles));  % Search for the same intersection in both directions
                if ~isempty(intersection) % If you have passed through the layer, this will not be empty
                    ix=triangles(intersection); % Recover ID of intersection triangle
                    ix=faces(ix,:); % Get corners of said triangle
                    pt(end,6)=mean([seis_depths(ix(1),3),seis_depths(ix(2),3),seis_depths(ix(3),3)]); % change to weighted mean
                    pt(end,7)=step;
                    trisurf(faces,x,y,z_down, intersectd*1.0,'FaceAlpha', 0.9);
%                     end_search=1;
                end
            end
        end
    end
    pt(step+1,1:3)=pt(step,1:3)+dir(step+1,:);
    if rem(step,5)==0
       plot3(pt(end,1),pt(end,2),pt(end,3),'rd','MarkerFaceColor','r')
    else
    plot3(pt(end,1),pt(end,2),pt(end,3),'r.')
    end
%     drawnow
    
    step=step+1;
    depth=pt(end,3);
    if pt(end,1)>250 || pt(end,2) < -100 || end_search==1
        break
    end
end

[P,temp]=seismic_PT(seis_depths,pt);
fprintf('Max Depth: %.1f km \n Fracture Time: %.2f Ma \n Max Temp: %.0f C\n Total Time: %.2f Ma\n',-pt(end,3),(pt(end,7)-pt(end,5))*timestep*1e-6,max(temp),step*timestep*1e-6)
%%

figure
plot3(pt(:,1),pt(:,2),pt(:,3),'.')
xlabel('X');ylabel('Y');zlabel('Z')
axis equal
grid
hold on
plot3(p([1:4 1],1),p([1:4 1],2),p([1:4 1],3))
plot3(p(1,1),p(1,2),p(1,3),'r.')
plot3(p(2,1),p(2,2),p(2,3),'b.')
plot3(p(3,1),p(3,2),p(3,3),'g.')
plot3(p(4,1),p(4,2),p(4,3),'k.')

for ii=1:step
    [m]=rotate3d(dir(ii,:)); % Make 3D rotation matrix from vector
    p(4*ii-3:4*ii,:)=(m\(p(1:4,:)-pt(1,1:3))')'+pt(ii,1:3); % Rotate original horizontal square, and translate to new center
    plot3(p([4*ii-3:4*ii 4*ii-3],1),p([4*ii-3:4*ii 4*ii-3],2),p([4*ii-3:4*ii 4*ii-3],3))
    plot3(p(4*ii-3,1),p(4*ii-3,2),p(4*ii-3,3),'r.')
    plot3(p(4*ii-2,1),p(4*ii-2,2),p(4*ii-2,3),'b.')
    plot3(p(4*ii-1,1),p(4*ii-1,2),p(4*ii-1,3),'g.')
    plot3(p(4*ii,1),p(4*ii,2),p(4*ii,3),'k.')
%     drawnow
end

%%
p(1:4,:)=p(5:8,:)-pt(2,1:3)+pt(1,1:3);
shp=alphaShape(p(:,1),p(:,2),p(:,3),4);
% shp=alphaShape(p(1:end-4,1),p(1:end-4,2),p(1:end-4,3),4);
a = criticalAlpha(shp,'one-region');
% a=ceil(shp.Alpha*10)/10;
figure
plot(shp)
axis equal
%%
gifpath(0,dir,timestep,pt,exhum,temp,P,gps,step,p,shp)