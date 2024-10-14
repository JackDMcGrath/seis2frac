

warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')
%% Prepare Velocity Fields
bearing=56.8962; % Fault Bearing
origin=[168.685,-44.115,0]; % UTM Origin
fault=[0,0;324.2663,0]; % UTM Fault coords
timestep=1e5; % timestep for path calculations
grid_res=1; % km

lon=-50:grid_res:250; % Define volume of interest
lat=-100:grid_res:20;
depth=-30:0;
[xx,yy,zz]=meshgrid(lon,lat,depth); % create coordinates

[vx,vy,vz]=velocity_fields(bearing,origin,xx,yy); % Fill this region with velocity fields using GPS and exhumation

vx=-vx;% vx(:)=0;
vy=-vy;
vz=-vz;

% Run back trace

% pt=[142.0609,-7.4839,0,0,0,0,0]; % Crawford's Knob
% pt=[141.6346,-2.6302,0,0,0,0,0]; % Franz Josef Planar

[sx,sy,sz]=meshgrid([140],[10],0); % Start points
s=stream3(xx,yy,zz,vx,vy,vz,sx,sy,sz,[1]); % Ray path

for ii=1:size(s,2) % Remove entries too close to the edge (no room for 1 step)
    ix=size(s,2)-ii+1;
    if size(s{ix},1)==1
        s(ix)=[];
    end
end
%%
f=figure('Name','With Fault Parallel')
h=streamline(s);
hold on
xlim([100 250])
ylim([-100 50])
pbaspect([1 1 0.5])
% plot(ll_rot(:,1),ll_rot(:,2),'k')
fill3(af([1 2 4 3 1],1),af([1 2 4 3 1],2),af([1 2 4 3 1],3),'r')
for ii=1:size(s,2)
    plot3(s{ii}(:,1),s{ii}(:,2),s{ii}(:,3),'.')
end
xlabel('X');ylabel('Y');zlabel('Depth');
view(3)



%% Prepare depth surfaces

load data/depth_surfaceTmean.mat;

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

poly_width=5; %width of poly in km

dir=s; % Make variable to hold vectors

X2=get(h,'XData'); %https://itectec.com/matlab/matlab-how-to-get-velocity-values-for-stream3-function/
Y2=get(h,'YData');
Z2=get(h,'ZData');

for ii=1:size(s,2)
    dir{ii}(1,:)=[0,0,0];
    clear p
    [p{ii}(1,:),p{ii}(2,:),p{ii}(3,:),p{ii}(4,:)]=find_corners(s{ii}(1,1:3),poly_width);
    for jj=2:size(s{ii},1)
        dir{ii}(jj,:)=s{ii}(jj,:)-s{ii}(jj-1,:); % Calculate vector
        [m]=rot3d(dir{ii}(jj,:)); % Make 3D rotation matrix from vector
        p{ii}(4*jj-3:4*jj,:)=(m\(p{ii}(1:4,:)-s{ii}(1,1:3))')'+s{ii}(jj,1:3); % Rotate original horizontal square, and translate to new center
    end
    
    if size(s,2)>1
        u1 = interp3(xx,yy,zz,vx,X2{ii},Y2{ii},Z2{ii});
        v1 = interp3(xx,yy,zz,vy,X2{ii},Y2{ii},Z2{ii});
        w1 = interp3(xx,yy,zz,vz,X2{ii},Y2{ii},Z2{ii});
    else
        u1 = interp3(xx,yy,zz,vx,X2,Y2,Z2);
        v1 = interp3(xx,yy,zz,vy,X2,Y2,Z2);
        w1 = interp3(xx,yy,zz,vz,X2,Y2,Z2);
    end
    s{ii}(1,4:6)=[0 0 0];
    s{ii}(:,4)=sqrt(u1.^2+v1.^2+w1.^2); % Speed
    s{ii}(2:end,5)=sqrt((s{ii}(2:end,1)-s{ii}(1:end-1,1)).^2+(s{ii}(2:end,2)-s{ii}(1:end-1,2)).^2+(s{ii}(2:end,3)-s{ii}(1:end-1,3)).^2); % Distance (km)
    s{ii}(:,6)=s{ii}(:,5)./s{ii}(:,4); % Time (yr) S=D/T -> T = D/S
end

fprintf('Max Depth: %.1f km \nTime: %.2f Ma\n',max(abs(s{1}(:,3))),sum(s{1}(:,6))*1e-6)
fprintf('Mean Velocity: %.1f mm/yr\n',(sqrt(sum((s{1}(end,1:3)-s{1}(1,1:3)).^2))*1e6)/sum(s{1}(:,6)))
%%
% depth=pt(end,3);
% max_depth=30;
%
% timestep=1e5;
% step=1;



% while depth<max_depth

%     [dx]=griddata(gps(:,1),gps(:,2),gps(:,3),pt(end,1),pt(end,2));
%     [dy]=griddata(gps(:,1),gps(:,2),gps(:,4),pt(end,1),pt(end,2));
%     [dz]=griddata(exhum(:,1),exhum(:,2),exhum(:,3),pt(end,1),pt(end,2));
%     [dx_slip]=griddata(xx(:),yy(:),sliprate_grid(:),pt(end,1),pt(end,2));
%     [dy_slip]=griddata(xx(:),yy(:),dipn_grid(:),pt(end,1),pt(end,2));
%     [dz_slip]=griddata(xx(:),yy(:),dipv_grid(:),pt(end,1),pt(end,2));
%     if isnan(dz)
%         dz=0.1;
%     end
%     dir(step+1,:)=-[dx*1e-6*timestep,dy*1e-6*timestep,dz*1e-6*timestep];
% %     dir(step+1,:)=-[dx*1e-6*timestep,dy*1e-6*timestep,(dz+dz_slip)*1e-6*timestep]; % Using the dipslip
%     dir(step+1,:)=-[(dx+dx_slip)*1e-6*timestep,(dy+dy_slip)*1e-6*timestep,dz*1e-6*timestep]; % Using the sliprates
%     dir(step+1,:)=-[(dx+dx_slip)*1e-6*timestep,(dy+dy_slip)*1e-6*timestep,(dz+dz_slip)*1e-6*timestep]; % Using the sliprates + dipslip

gps=readmatrix('data/gps.csv');
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

exhum=readmatrix('C:\Jacko\NZ 2020\seis2frac3\data\Modeled_exhumation_rates_grid_points.csv');
exhum([11,1:10],3)=exhum(12:22,3); %Make AUS nodes equal PACmax
exhum(:,(1:2))=ll2utm(exhum(:,(1:2)),origin(:,1:2));
[exrot]=point_rotation(exhum(:,1:2),90-bearing,'strike');
exhum(:,1:2)=exrot(:,1:2);

for ii=1:size(s,2)
    figure;
    hold on
    xlabel('X');ylabel('Y');zlabel('Z')
    % plot3(pt(end,1),pt(end,2),pt(end,3),'r.');
    view([-45,30])
    grid
    
    searchupper=1;
    searchlower=0;
    end_search=0;
    intersectusurf=zeros(size(vu1,1),1);
    intersectdsurf=zeros(size(vu1,1),1);
    for jj=1:size(s{ii},1)-1
        
        if searchupper==1
            intersectu = TriangleRayIntersection(s{ii}(jj,1:3), dir{ii}(jj+1,1:3), vu1, vu2, vu3,'LineType','Segment');
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
                %                     searchupper=0;
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
            intersectd = TriangleRayIntersection(s{ii}(jj,1:3), dir{ii}(jj+1,1:3), vd1, vd2, vd3,'LineType','segment');
            if ~isempty(find(intersectd))
                %             backcheck = TriangleRayIntersection(pt(end,1:3)+dir(step+1,:), -dir(step+1,:), vd1, vd2, vd3);
                %             if ~isempty(find(backcheck)) % If the return vector also passes the surface
                %                 triangles=find(intersectd); % Find intersection of forward vector
                %                 back_triangles=find(backcheck); % Find intersetion of back vector
                %                 intersection=find(ismember(triangles,back_triangles));  % Search for the same intersection in both directions
                %                 if ~isempty(intersection) % If you have passed through the layer, this will not be empty
                %                     ix=triangles(intersection); % Recover ID of intersection triangle
                %                     ix=faces(ix,:); % Get corners of said triangle
                %                     pt(end,6)=mean([seis_depths(ix(1),3),seis_depths(ix(2),3),seis_depths(ix(3),3)]); % change to weighted mean
                %                     pt(end,7)=step;
                intersectdsurf(find(intersectd))=1;
                trisurf(faces,x,y,z_down, intersectdsurf*1.0,'FaceAlpha', 0.9);
                %                     end_search=1;
                %                 end
                %             end
            end
        end
        %     pt(step+1,1:3)=pt(step,1:3)+dir(step+1,:);
        %     if rem(step,5)==0
        %        plot3(pt(end,1),pt(end,2),pt(end,3),'rd','MarkerFaceColor','r')
        %     else
        %     plot3(pt(end,1),pt(end,2),pt(end,3),'r.')
        %     end
        %     drawnow
        
        %     step=step+1;
        %     depth=pt(end,3);
        %     if pt(end,1)>250 || pt(end,2) < -100 || end_search==1
        %         break
        %     end
        
        plot3(s{ii}(:,1),s{ii}(:,2),s{ii}(:,3),'r.')
    end
    
    [P,temp]=seismic_PT(seis_depths,s{ii});
    fprintf('Max Depth: %.1f km \n Fracture Time: %.2f Ma \n Max Temp: %.0f C\n',max(abs(s{ii}(:,3))),0,max(temp))
    
    
    %
    
    poly_width=5; %width of poly in km
    p=nan(4,3);
    [p(1,:),p(2,:),p(3,:),p(4,:)]=find_corners(s{ii}(1,1:3),poly_width);
    
    % figure
    % plot3(pt(:,1),pt(:,2),pt(:,3),'.')
    % xlabel('X');ylabel('Y');zlabel('Z')
    % axis equal
    % grid
    % hold on
    % plot3(p([1:4 1],1),p([1:4 1],2),p([1:4 1],3))
    % plot3(p(1,1),p(1,2),p(1,3),'r.')
    % plot3(p(2,1),p(2,2),p(2,3),'b.')
    % plot3(p(3,1),p(3,2),p(3,3),'g.')
    % plot3(p(4,1),p(4,2),p(4,3),'k.')
    
    
    for jj=1:size(s{ii},1)
        [m]=rot3d(dir{ii}(jj,:)); % Make 3D rotation matrix from vector
        p(4*jj-3:4*jj,:)=(m\(p(1:4,:)-s{ii}(1,1:3))')'+s{ii}(jj,1:3); % Rotate original horizontal square, and translate to new center
        %     plot3(p([4*ii-3:4*ii 4*ii-3],1),p([4*ii-3:4*ii 4*ii-3],2),p([4*ii-3:4*ii 4*ii-3],3))
        %     plot3(p(4*ii-3,1),p(4*ii-3,2),p(4*ii-3,3),'r.')
        %     plot3(p(4*ii-2,1),p(4*ii-2,2),p(4*ii-2,3),'b.')
        %     plot3(p(4*ii-1,1),p(4*ii-1,2),p(4*ii-1,3),'g.')
        %     plot3(p(4*ii,1),p(4*ii,2),p(4*ii,3),'k.')
        %     drawnow
    end
    
    %%
    
    
    p(1:4,:)=p(5:8,:)-s{ii}(2,1:3)+s{ii}(1,1:3);
    shp=alphaShape(p(:,1),p(:,2),p(:,3),4);
    a = criticalAlpha(shp,'one-region');
    % % a=ceil(shp.Alpha*10)/10;
    % figure
    % plot(shp)
    % axis equal
    %
    
    try
        gifpath(0,dir{ii},timestep,s{ii},exhum,temp,P,gps,size(s{ii},1),p,shp)
    catch
        close
    end
end