function [ray,path_poly,ftime,pressure,temp,vels,AF]=exhumation_path(stepsize,grid_res,fault_base,origin,bearing,raystart,gpsfile,exhum,seis_depths,af,poly_width,forward_ray,compress)

warning('off','MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId')

% Prepare Velocity Fields

lon=0:grid_res:200; % Define volume of interest
lat=-50:grid_res:0;
depth=0:-1:-fault_base;
[xx,yy,zz]=meshgrid(lon,lat,depth); % create coordinates

%[vx,vy,vz]=velocity_fields(bearing,origin,xx,yy,gpsfile,exhum); % Fill this region with velocity fields using GPS and exhumation

strikeSlip=-27e-6*0; % km/yr
converge=10e-6; % km/yr
dipSlip=8e-6; % km/yr, dip slip on the fault
[vx,vy,vz,AF]=velocity_fields_new(bearing,origin,xx,yy,zz,strikeSlip,converge,dipSlip,exhum,af); % Fill this region with velocity fields using GPS and exhumation
afSurface=AF.surface;
velTot=sqrt(vx.^2+vy.^2+vz.^2); % Velocity magnitude in volume

% Reverse the sense of the velocity when backprojecting rays
if forward_ray==0
    vx=-vx;% vx(:)=0;
    vy=-vy;
    vz=-vz;
end

% Run back trace iteritively from base up, allowing for rock interaction

% [sx,sy,sz]=meshgrid(raystart(:,1),raystart(:,2),raystart(:,3)); % Start points MAY BE ABLE TO ACCOMODATE MULTIPLE?
% ray=stream3(xx,yy,zz,vx,vy,vz,sx,sy,sz,stepsize); % Ray path

depths=unique(raystart(:,3));

ray={};

% Create intersection surface for AF
AF.faces=boundary(AF.surface(:,1),AF.surface(:,2),AF.surface(:,3));
% Verticies of each delaunay triangulation
af1=afSurface(AF.faces(:,1),:);
af2=afSurface(AF.faces(:,2),:);
af3=afSurface(AF.faces(:,3),:);

if numel(unique(raystart(:,1)))==1
    width=10;
else
    xstarts=unique(raystart(:,1));
    width=abs(xstarts(2)-xstarts(1));
end

if numel(depths)==1
    height=10;
else
    height=abs(depths(2)-depths(1));
end


plot_figure=0;
if plot_figure==1
figure
plot3(AF.surface(:,1),AF.surface(:,2),AF.surface(:,3),'r*')
hold on
plot3(AF.surface(:,1),AF.surface(:,2),AF.surface(:,3))
view([-90 0])
end

for ii=1:numel(depths)
    ix=find(raystart(:,3)==depths(ii));
    newray=stream3(xx,yy,zz,vx,vy,vz,raystart(ix,1),raystart(ix,2),raystart(ix,3),stepsize); % Ray path
    for jj=1:numel(ix)
        if ~isempty(newray{jj})
            dir{jj}=newray{jj}(:,1:3); % Make variable to hold vectors
            dir{jj}(1,:)=[0,0,0];
            [path_poly{jj}(1,:),path_poly{jj}(2,:),path_poly{jj}(3,:),path_poly{jj}(4,:)]=find_corners(newray{jj}(1,1:3),poly_width);
            for kk=2:size(newray{jj},1)
                dir{jj}(kk,:)=newray{jj}(kk,:)-newray{jj}(kk-1,:); % Calculate vector
            end
            
            % Search each ray and see if one has intersected AF region
            if plot_figure==1
            plot3(newray{jj}(:,1),newray{jj}(:,2),newray{jj}(:,3),'k*')
            end
            for kk=1:size(newray{jj},1)-1
                intersectAF = TriangleRayIntersection(newray{jj}(kk,1:3), dir{jj}(kk+1,1:3), af1, af2, af3,'LineType','Segment');
                if plot_figure==1
                plot3(newray{jj}(kk,1),newray{jj}(kk,2),newray{jj}(kk,3),'g*')
                plot3(newray{jj}(kk+1,1),newray{jj}(kk+1,2),newray{jj}(kk+1,3),'b*')
                end
                if ~isempty(find(intersectAF, 1))
                    % Ray has entered AF velocity region
                    fprintf('Depth %.0f km Ray %.0f has entered AF region at step %.0f\nAdjusting velocity field %.1f km\n',depths(ii),jj,kk, height/compress)
                    [vx,vy,vz,AF]=updateVelocityField(xx,yy,zz,vx,vy,vz,AF,poly_width,height,compress,newray{jj}(kk+1,:)');
                    % Create intersection surface for AF
                    % Verticies of each delaunay triangulation
                    af1=AF.surface(AF.faces(:,1),:);
                    af2=AF.surface(AF.faces(:,2),:);
                    af3=AF.surface(AF.faces(:,3),:);
                    break
                end
            end
        end
        
        ray{size(ray,2)+1}=newray{jj};
        
    end
    
end


for ii=1:size(ray,2) % Remove entries too close to the edge (no room for 1 step)
    ix=size(ray,2)-ii+1;
    if size(ray{ix},2)==1
        ray(ix)=[];
    end
end

%% Prepare depth surfaces

x=seis_depths(~isnan(seis_depths(:,1)),1);
y=seis_depths(~isnan(seis_depths(:,1)),2);
faces=delaunay(x,y); % Triangulate seis_depths

z_up=seis_depths(~isnan(seis_depths(:,1)),3); % Top surface depths
z_down=seis_depths(~isnan(seis_depths(:,1)),4); % Lower Surface Depths

v_up=[x,y,z_up]; % Top Surface 3D co-ordinates
v_down=[x,y,z_down]; % Lower surface 3D co-ordinates

% Verticies of each delaunay triangulation
vu1=v_up(faces(:,1),:);
vu2=v_up(faces(:,2),:);
vu3=v_up(faces(:,3),:);

vd1=v_down(faces(:,1),:);
vd2=v_down(faces(:,2),:);
vd3=v_down(faces(:,3),:);

path_poly=cell(1,size(ray,2)); % Variable to hold the polygon points around the the raypath % NO LONGER NEEDED IF USING DATADENSITY

for ii=1:size(ray,2)
    % Start Firing Rays through the seismicity levels
    if ~isempty(ray{ii})
        dir{ii}=ray{ii}(:,1:3); % Make variable to hold vectors
        dir{ii}(1,:)=[0,0,0];
        [path_poly{ii}(1,:),path_poly{ii}(2,:),path_poly{ii}(3,:),path_poly{ii}(4,:)]=find_corners(ray{ii}(1,1:3),poly_width);
        for jj=2:size(ray{ii},1)
            dir{ii}(jj,:)=ray{ii}(jj,:)-ray{ii}(jj-1,:); % Calculate vector
            [m]=rot3d(dir{ii}(jj,:)); % Make 3D rotation matrix from vector
            path_poly{ii}(4*jj-3:4*jj,:)=(m\(path_poly{ii}(1:4,:)-ray{ii}(1,1:3))')'+ray{ii}(jj,1:3); % Rotate original horizontal square, and translate to new center
        end
        
        ray{ii}(2:size(ray{ii},1),4)=interp3(xx,yy,zz,velTot,ray{ii}(1:end-1,1),ray{ii}(1:end-1,2),ray{ii}(1:end-1,3)); % Vel to reach each point
        ray{ii}(:,5)=[0;sqrt(sum((ray{ii}(1:end-1,:)-ray{ii}(2:end,:))'.^2)')]; % Distance between each step
        
        ray{ii}(:,6)=[0;ray{ii}(2:end,5)./ray{ii}(2:end,4)]; % Timestep
        
        plot_figs=0;
        if plot_figs==1
            figure;
            hold on
            xlabel('Along Strike (km)');ylabel('Cross Strike (km)');zlabel('Depth (km)')
            title('Upper and Lower Seismicity Surfaces')
            trisurf(faces,x,y,z_up);
            trisurf(faces,x,y,z_down);
            fill3(af([1 2 4 3 1],1),af([1 2 4 3 1],2),af([1 2 4 3 1],3),'g')
            view([-45,30])
            grid
        end
        
        if forward_ray==0
            searchupper=1;
            searchlower=0;
            intersectusurf=zeros(size(vu1,1),1);
            intersectdsurf=zeros(size(vu1,1),1);
            for jj=1:size(ray{ii},1)-1
                if plot_figs==1
                    plot3(ray{ii}(jj,1),ray{ii}(jj,2),ray{ii}(jj,3),'r.')
                end
                if searchupper==1
                    intersectu = TriangleRayIntersection(ray{ii}(jj,1:3), dir{ii}(jj+1,1:3), vu1, vu2, vu3,'LineType','Segment');
                    if ~isempty(find(intersectu, 1))
                        searchlower=1;
                        intersectusurf(intersectu==1)=1;
                        if plot_figs==1
                            trisurf(faces,x,y,z_up, intersectusurf*1.0,'FaceAlpha', 0.9);
                        end
                    end
                end
                
                if searchlower==1
                    intersectd = TriangleRayIntersection(ray{ii}(jj,1:3), dir{ii}(jj+1,1:3), vd1, vd2, vd3,'LineType','segment');
                    if ~isempty(find(intersectd, 1))
                        intersectdsurf(intersectd==1)=1;
                        if plot_figs==1
                            trisurf(faces,x,y,z_down, intersectdsurf*1.0,'FaceAlpha', 0.9);
                        end
                    end
                end
            end
        else
            searchupper=0;
            searchlower=1;
            intersectusurf=zeros(size(vu1,1),1);
            intersectdsurf=zeros(size(vu1,1),1);
            for jj=1:size(ray{ii},1)-1
                if plot_figs==1
                    plot3(ray{ii}(jj,1),ray{ii}(jj,2),ray{ii}(jj,3),'r.')
                end
                if searchupper==1
                    intersectu = TriangleRayIntersection(ray{ii}(jj,1:3), dir{ii}(jj+1,1:3), vu1, vu2, vu3,'LineType','Segment');
                    if ~isempty(find(intersectu, 1))
                        intersectusurf(intersectu==1)=1;
                        if plot_figs==1
                            trisurf(faces,x,y,z_up, intersectusurf*1.0,'FaceAlpha', 0.9);
                        end
                    end
                end
                
                if searchlower==1
                    intersectd = TriangleRayIntersection(ray{ii}(jj,1:3), dir{ii}(jj+1,1:3), vd1, vd2, vd3,'LineType','segment');
                    if ~isempty(find(intersectd, 1))
                        intersectdsurf(intersectd==1)=1;
                        searchupper=1;
                        if plot_figs==1
                            trisurf(faces,x,y,z_down, intersectdsurf*1.0,'FaceAlpha', 0.9);
                        end
                    end
                end
            end
        end
        
        % Calculate pressures, temperatures and fracture time
        [pressure{ii},temp{ii},ftime{ii}]=seismic_PT(seis_depths,ray{ii});
        %%
        fprintf('\n*** Ray %.0f ***\n Mean Velocity: %.1f mm/yr\n Max Depth: %.1f km \n Total Time: %.2f Ma\n Fracture Time: %.2f Ma \n Max Temp: %.0f C\n Timesteps: %.0f\n', ...
            ii,(sqrt(sum((ray{ii}(end,1:3)-ray{ii}(1,1:3)).^2))*1e6)/sum(ray{ii}(:,6)),max(abs(ray{ii}(:,3))),sum(ray{ii}(:,6))*1e-6,nansum(ray{ii}(ftime{ii},6))/1e6,max(temp{ii}),size(ftime{ii},1))
        
    end
    
vels.vx=vx;
vels.vy=vy;
vels.vz=vz;
vels.xx=xx;
vels.yy=yy;
vels.zz=zz;
    
end
