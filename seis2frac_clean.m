%% SEIS2FRAC
% Script for analysing seismicity and relating to fracture propeties -using
% just quadtree subsampling and backtracing exhumation path
% https://tinyurl.com/3s7tmyfe

%% Setup
clear
addpath(genpath('seis2frac'))
addpath(genpath('C:/Jacko/scripts'))
% Toggles
plot_figures=0;
forward_ray=1; % Toggle for if we start at the surface and work backwards, or vice versa for the ray paths. Working forwards intended to allow us to interact with the fault

% File paths
eq_files={'data\seismic_data.csv';'data\DWARFS_North_Catalog_Final_Revised.csv'}; %'data\DWARFS_South_Catalog_Final_Revised.csv'};
coast_file='C:\Jacko\scripts\coastlines\gshhg-bin-2.3.7';
exhum_file='C:\Jacko\NZ 2020\seis2frac3\data\Modeled_exhumation_rates_grid_points.csv';

% Quadtree Geometry
gridx=10; % set size of cube or minSize of Quadtree (km)
QuadminbinCapacity=30;
quadtree_iterations=2; % Number of times to run quadtree algorythm to deal with unbinned data

% Quadtree Rotation flags
r_flg.norm2strike=1; % Rotate dataset to strike before quadtree
r_flg.norm2dip=0; % Rotate dataset to dip before quadtree
r_flg.strike2norm=1; % Rotate dataset back from strike to normal after quadtree
r_flg.dip2norm=0; % Rotate dataset back from dip to normal after quadtree

% Fault Geometry
af_ll=[168.685,-44.115,0;171.9880,-42.4733,0]; % Lon Lat of Fault Trace
dip=45; % set fault dip (degrees)
fault_base=35; % set depth of fault (km)

% Rock Properties
mu=22e9; % Shear modulus, Pa
vein_thickness=1.5; % mm
vein_lengths=[1 10]; % m

% Fault Length-Displacement Scaling Parameters
scaling_a=0.0337; % 3e-5
scaling_exp=1.02; % 1.0

% G-R Constants from Michalios 2019
avalue=[]; % If empty, will calculate from dataset (4.28)

if size(eq_files,1)==1 && strcmpi('data\seismic_data.csv',eq_files{1})
    bvalue=0.85; % If empty, will calculate from dataset. Value from 9111 events =0.85 (Michalios et al.)
else
    bvalue=[];
end

MOC=1.1; % Magnitude of Completeness

%% Prepare Data
% Fault Limits
af=af_ll; % Variable to change to UTM

if ~isempty(af)
    origin=af(1,1:2); % Set UTM origin as one end of fault
    af(:,1:2)=ll2utm(af(:,1:2),origin); % Set to UTM
    bearing=atand((af(2,1)-af(1,1))/(af(2,2)-af(1,2))); % Work out fault bearing
    
    % Calculation: Co-ordinates of Fault Base
    af(3:4,:)=[af(1,1)+(cosd(bearing)*(fault_base/tand(dip))),af(1,2)-(sind(bearing)*(fault_base/tand(dip))),-fault_base;...
        af(2,1)+(cosd(bearing)*(fault_base/tand(dip))),af(2,2)-(sind(bearing)*(fault_base/tand(dip))),-fault_base];
else
    origin=mean(eqs(:,[1 2]));
end

% Load Seismicity data
eqs=[];
eq_source=1;
% Some seismicity has values > 0km (ie in topography. Is this an issue?)
for ii=1:size(eq_files,1)
    eqs_tmp=readmatrix(eq_files{ii}); % Import seismicity data
    if size(eqs_tmp,2)==12
        eqs=[eqs;[eqs_tmp(:,8),eqs_tmp(:,7),-eqs_tmp(:,9),eqs_tmp(:,10)]]; % Lon, Lat, Depth (-ve), Mag
    else
        eqs=[eqs;[eqs_tmp(:,3),eqs_tmp(:,2),-eqs_tmp(:,4),eqs_tmp(:,6)]]; % Lon, Lat, Depth (-ve), Mag
    end
    eq_source=[eq_source,eq_source(end)+size(eqs_tmp,1)]; % Indicies of the start and end of new catalogue
end
clear eqs_tmp
eqs(:,5:6)=nan; % Lon, Lat, Depth, Mag, Dist to Fault, Moment
eqs(eqs(:,3)<-fault_base,:)=[]; % Remove any seismicity > fault base
eqs(:,(1:2))=ll2utm(eqs(:,(1:2)),origin); % Set to UTM
eqs(:,5)=point_to_line(eqs(:,(1:2)),af(1,:),af(2,:)); % Calculate horizontal distance to fault trace
eqs(:,6)=10.^((3/2).*(eqs(:,4)+6.07)); % Calculate Seismic Moment of each event
eqs(:,7)=0; % Label each eq with it's catalogue
for ii=1:size(eq_files,1)
    eqs(eq_source(ii):eq_source(ii+1)-1,7)=ii;
end

% Load Coastline data if needed
if ~exist('ll','var') && plot_figures==1
    [ll(:,1),ll(:,2)]=get_coast(168,173.25,-45,-41.75,[],'f',coast_file);
    ll_utm=ll2utm(ll(:,1:2),origin);
elseif ~exist('ll','var') && plot_figures==0
    ll=[];
    ll_utm=[];
end

%% Gutenburg-Richter Plots

% Set Histogram Parameters (0.1 Mw Bins)
bins=floor(min(eqs(:,4))*10)/10:0.1:ceil(max(eqs(:,4))*10)/10;

[n,~,~,btrend,GRmdl]=gutenberg_richter(eqs(:,4),bins,bvalue,[],MOC,plot_figures);

if isempty(avalue)
    avalue=GRmdl.Coefficients.Estimate(1);
end

if isempty(bvalue)
    bvalue=-GRmdl.Coefficients.Estimate(2);
end

%% Calculate displacement of fractures MAY BE WORTH SPLITTING INTO SEPERATE FUNCTION SOON ENOUGH
% https://www.sciencedirect.com/science/article/pii/S0191814119303657
% dmax=0.0337L^1.02
% I am starting to suspect the validity of this law (specifically the use
% of Dmax). maybe try disp=3e-5 L (Rich's paper
% https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/jgrb.50236 or https://www.sciencedirect.com/science/article/pii/S0040195113003776)

[dis,mag_window,len]=fracmag2disp(bins,scaling_a,scaling_exp,mu,vein_lengths);

fprintf('For Fractures %.1f - %.1f m, magnitudes = Mw %.1f to %.1f\n', vein_lengths(1),vein_lengths(2),mag_window(1),mag_window(2))

%% Quadtree subsampling of the data

dropped=cell(1,quadtree_iterations);

[QT,dropped{1}] = quadtree_subsample(eqs(:,1:4),af,bearing,dip,QuadminbinCapacity,r_flg,ll_utm,'shrink',plot_figures,'centeredSquare',50);

%% Fill gaps by iterating twice more onto unbinned seismicity

if quadtree_iterations>1
    eqd=eqs; % seismicity variable that will be cut down iteratively
    
    for ii=1:quadtree_iterations-1
        
        [DT,dropped{ii+1},llrot] = quadtree_subsample(eqd(dropped{ii},1:4),af,bearing,dip,QuadminbinCapacity/2,r_flg,ll_utm,'shrink',plot_figures,'centeredSquare',50);
        
        % Add newly QuadTree-d data to main QT
        QT.Points=[QT.Points;DT.Points];
        QT.PointBins=[QT.PointBins;DT.PointBins+QT.BinCount];
        QT.BinCount=QT.BinCount+DT.BinCount;
        QT.BinBoundaries=[QT.BinBoundaries;DT.BinBoundaries];
        QT.BinCorners=cat(3,QT.BinCorners,DT.BinCorners);
        
        eqd=eqd(dropped{ii},1:3);
    end
    
    % Find any seismicity that still hasn't been binned
    for ii=quadtree_iterations-1:-1:1
        ix=(1:length(dropped{ii}))';
        ix=ismember(ix,dropped{ii+1});
        dropped{ii}(ix==0)=[];
    end
    
    eqs(dropped{1},:)=[]; % Remove unbinned data from eqs variable
    
    clear eqd DT dropped ix ii
end

if plot_figures==1
    figure
    figname='QuadTree Subsample - Final';
    
    set(gcf,'renderer','zbuffer','name',figname); title(figname);
    plot(ll(:,1),ll(:,2))
    hold on
    boxH = QT.plot;
    cols = lines(QT.BinCount);
    doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
    for i = 1:QT.BinCount
        set(boxH(i),'Color',cols(i,:),'LineWidth', 1)
        doplot3(QT.Points(QT.PointBins==i,:),'.','Color',cols(i,:))
    end
    xlabel('Lon');ylabel('Lat');zlabel('Depth');
    plot3(117.4818,56.1926,3.5,'kp','MarkerFaceColor','g','MarkerEdgeColor','w','MarkerSize',12)
    plot(llrot(:,1),llrot(:,2))
    plot3(af([1 2 4 3 1],1),af([1 2 4 3 1],2),af([1 2 4 3 1],3),'LineWidth',3,'Color','k');
    xlim([-50 300]);ylim([-100 200]);
    hold off
    axis equal
    xlim([floor(min(QT.Points(:,1))/50)*50 ceil(max(QT.Points(:,1))/50)*50]);ylim([floor(min(QT.Points(:,2))/50)*50 ceil(max(QT.Points(:,2))/50)*50])
end

% Problem with this is that some of the drop iterations cover the same area
% as other bins. May need to assign a heirarchy to these bins in the future


%% Determine Upper and Lower Limits of Seismicity
% NOTE: ALTERING THIS SECTION CURRENTLY DOES NOTHING AS THE NEXT SECTION
% WILL CALL A VARIABLE THAT IS A GMT SURFACE PREVIOUSLY CALCULATED FROM
% THIS POINTS. NO I DON'T KNOW HOW I DID THAT

seis_depths=nan(QT.BinCount,4);
bin_center=nan(QT.BinCount,2);

for ii=1:QT.BinCount
    seis_depths(ii,1:2)=mean(QT.BinCorners(1:4,1:2,ii)); % Find bin center
    if sum(QT.PointBins==ii)>=QuadminbinCapacity
        [seis_depths(ii,3),seis_depths(ii,4)] = upper_lower_seis(QT.Points(find(QT.PointBins==ii),3),0);
    end
    bin_center(ii,:)=mean(QT.BinCorners(1:4,1:2,ii));
end

seis_ix=find(~isnan(seis_depths(:,1))); % identify bins that have upper-lower depths


%% Exhumation path time

% Load exhumation data
exhum=readmatrix(exhum_file);
exhum(:,(1:2))=ll2utm(exhum(:,(1:2)),origin);

exhum(1:11,3)=exhum(12:22,3);
fprintf('NB: AUS Plate nodes copied to max\n')


grid_res=1; % km
stepsize=1; % path calculation step size(in terms of grid_res. Note- this is n-steps in x,y,or z, not vector magnitude)

alongStrikeStarts=(100:10:150);
crossStrikeStarts=(-50:10:0 );
depthRange=(-35:5:-5);

gpsfile='data/gps.csv';
gmtsurfaces=load('data/depth_surfaceTmean.mat');
gmtsurfaces=gmtsurfaces.seis_depths;

[af_rot(:,1:3)]=point_rotation(af(:,1:2),90-bearing,'strike');
af_rot(:,3)=af(:,3);

[eqs_rot(:,1:3)]=point_rotation(eqs(:,1:2),90-bearing,'strike');
eqs(:,1:2)=eqs_rot(:,1:2);
eqs(eqs(:,4)<mag_window(1),:)=[];
eqs(eqs(:,4)>mag_window(2),:)=[];

poly_width=10;
volume_change=1/0.7;
volume_change=2;

if forward_ray==0
    [XX,YY]=meshgrid(alongStrikeStarts,crossStrikeStarts);
    raystart=[XX(:),YY(:)];
    raystart(:,3)=0;
    
else
    start_dist=-50; % Fault perpendicular distance from the fault to start
    [XX,ZZ]=meshgrid(alongStrikeStarts,depthRange);
    raystart=zeros(numel(XX),3)+start_dist;
    raystart(:,[1,3])=[XX(:),ZZ(:)];
end
%%
[ray,pathpoly,ftime,PT.pressure,PT.temp,vels,AF]=exhumation_path(stepsize,grid_res,fault_base,origin,bearing,raystart,gpsfile,exhum,gmtsurfaces,af_rot,poly_width,forward_ray,volume_change);
%%
figure('Name','Ray Paths')
hold on
for ii=1:size(ray,2)
    if size(ray{ii},2)>1
        plot3(ray{ii}(:,1),ray{ii}(:,2),ray{ii}(:,3))
    end
end
for ii=1:size(eq_files,1)
    plot3(eqs(eqs(:,7)==ii,1),eqs(eqs(:,7)==ii,2),eqs(eqs(:,7)==ii,3),'.')
end
fill3(af_rot([1 2 4 3 1],1),af_rot([1 2 4 3 1],2),af_rot([1 2 4 3 1],3),'g')
axis equal
view([-25 25])

%%
transect_density=zeros(size(raystart,1),1);
disvol_total=ones(size(raystart,1),1);
fracture_total_volume=zeros(size(raystart,1),1);
fracture_ratio=zeros(size(raystart,1),1);
eq_tot=nan(size(raystart,1),1);
eqMOC=(eqs(:,4)>=MOC); % Indicies of all events greater than MOC
observation_time=5; % Number of years that observations have been made

clear shp eqsin bin_btrend rock_volume eqs_total total_events fracture_volumes

method=2;

if method==1
    t='Individual Event Search';
    inpoly_transect_estimation % Method searches for individual events in raypath polygon
elseif method==2
    t='Data Density Search';
    t=sprintf('Surface Fracture Density (Compression = %.0f%%)',100-100/volume_change);
    data_density_transect_estimation % Method searches events volume polygon
elseif method==3
    t='Data Density w/ Variable Dissolution';
    variable_disolution_dd_transect_estimation % Method searches events volume polygon, and using variable amounts of dissolution
end
%% Plotting

load('Data/bamako.mat');
tmp=sort(transect_density(~isinf(transect_density(transect_density~=0))));
maxd=ceil(tmp(round(0.95*numel(tmp)))/5)*5;
try
    mind=min([0 floor(tmp(round(0.05*numel(tmp)))/5)*5]);
catch
    mind=0;
end


% Load Coastline data if needed
if ~exist('ll','var') || isempty(ll)
    [ll(:,1),ll(:,2)]=get_coast(168,173.25,-45,-41.75,[],'f',coast_file);
    ll_utm=ll2utm(ll(:,1:2),origin);
end

[ll_rot(:,1:3)]=point_rotation(ll_utm(:,1:2),90-bearing,'strike');

if forward_ray==0
    plotray=raystart;
    Xloc=XX;
    Yloc=YY;
else
    for ii=1:size(ray,2)
        if size(ray{ii},2)>1
        plotray(ii,:)=ray{ii}(end,1:3);
        end
    end
    
    Xloc=plotray(:,1);
    Yloc=plotray(:,2);
    
end

    
    
    scatter_plot(plotray,transect_density,[mind maxd],t,'Fractures/m',af_rot,ll_rot,bamako)
    ray_plot(ray,raystart,eqs,MOC,af_rot,ll_rot,dmap,dmapX,dmapY,dmapZ)
    title(sprintf('Ray Paths (Compression = %.0f%%)',100-100/volume_change))
    
    try
    contour_plot(plotray,transect_density,Xloc,Yloc,[mind maxd],1,25,t,'Fractures/m',af_rot,ll_rot,bamako);
    for ii=1:size(PT.temp,2)
        temp(ii)=max(PT.temp{ii});
        pressure(ii)=max(PT.pressure{ii});
        depth(ii)=-min(ray{ii}(:,3));
    end
    
    [contour_temp]=contour_plot(plotray,temp,Xloc,Yloc,[0 round(max(temp)/5)*5],50,100,sprintf('Peak Ray Temperatures (Compression = %.0f%%)',100-100/volume_change),'Temperature (C)',af_rot,ll_rot,bamako);
    %[contour_pressure]=contour_plot(plotray,pressure,Xloc,Yloc,[0 round(max(pressure)/5)*5],100,200,sprintf('Peak Ray Pressure (Compression = %.0f%%)',100-100/volume_change),'Pressure (MPa)',af_rot,ll_rot,bamako);
    [contour_depth]=contour_plot(plotray,depth,Xloc,Yloc,[0 ceil(max(depth)/5)*5],5,10,sprintf('Peak Depth (Compression = %.0f%%)',100-100/volume_change),'Depth (km)',af_rot,ll_rot,bamako);
    end
    
    figure
imagesc(vels.xx(1,:),vels.yy(:,1),squeeze(vels.vz(:,:,1)*1e6))
c=colorbar;
c.Label.String='mm/yr';
title(sprintf('Surface VZ rates (Compression = %.0f%%)',100-100/volume_change))
axis xy
caxis([0 8])

figure
trisurf(AF.faces,AF.surface(:,1),AF.surface(:,2),AF.surface(:,3))
title(sprintf('AF Velocity Regions (Compression = %.0f%%)',100-100/volume_change))
xlabel('Along Strike Distance (km)')
% ylabel('Cross Strike Distance (km)')
view([25 20])
zlabel('Depth (km)')