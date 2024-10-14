%% SEIS2FRAC
% Script for analysing seismicity and relating to fracture propeties -using
% just quadtree subsampling and backtracing exhumation path
% https://tinyurl.com/3s7tmyfe

%% Setup
% Admin
plot_figures=1;
eq_file='seismic_data.csv';
coast_file='C:\Jacko\scripts\coastlines\gshhg-bin-2.3.7';

% Geometry
gridx=10; % set size of cube or minSize of OcTree
QuadminbinCapacity=30;
dip=60; % set fault dip
fault_base=25; % set depth of fault

%Rock properties
mu=22e9; % Shear modulus, Pa
vein_thickness=1.5; % mm
vein_lengths=[1 10]; % m

% G-R Constants from Michalios 2019
avalue=[]; % If empty, will calculate from dataset (4.28)
bvalue=0.85; % If empty, will calculate from dataset. Value from 9111 events =0.85
MOC=1.1;

% Fault Length-Displacement Scaling Parameters
scaling_a=0.0337;
scaling_exp=1.02;
% scaling_a=3e-5;
% scaling_exp=1.0;

%% Prepare Data
% Seismic Catalogue
eqs=readmatrix(eq_file); %Import seismicity data
eqs=[eqs(:,8),eqs(:,7),-eqs(:,9),eqs(:,10)]; % Lon, Lat, Depth (-ve), Mag
eqs(:,[5:6])=nan; % Lon, Lat, Depth, Mag, Dist to Fault, Moment
eqs(eqs(:,3)<-fault_base,:)=[]; % Remove any seismicity > fault base
% Fault Limits
af_ll=[168.685,-44.115,0;171.9880,-42.4733,0]; % Variable to be kept in lon-lat
af=af_ll; % Variable to change to UTM
if ~isempty(af)
    origin=af(1,1:2); % Set UTM origin as one end of fault
    % Calculation: co-ords of base of fault
    af(:,1:2)=ll2utm(af(:,1:2),origin); % Set to UTM
    bearing=atand((af(2,1)-af(1,1))/(af(2,2)-af(1,2))); % Work out fault bearing
    % Calculation: Horizontal distance to Fault trace (Map view) - Lon Lart of fault base
    af(3:4,:)=[af(1,1)+(cosd(bearing)*(fault_base/tand(dip))),af(1,2)-(sind(bearing)*(fault_base/tand(dip))),-fault_base;...
        af(2,1)+(cosd(bearing)*(fault_base/tand(dip))),af(2,2)-(sind(bearing)*(fault_base/tand(dip))),-fault_base];
else
    origin=mean(eqs(:,[1 2]));
end

eqs(:,(1:2))=ll2utm(eqs(:,(1:2)),origin); % Set to UTM
% Calculation: Earthquake distance to fault
eqs(:,5)=point_to_line(eqs(:,(1:2)),af(1,:),af(2,:));
% Calculation: Seismic Moment of each event
eqs(:,6)=10.^((3/2).*(eqs(:,4)+6.07));

if ~exist('ll') && plot_figures==1
    [ll(:,1),ll(:,2)]=get_coast(168,173.25,-45,-41.75,[],'f',coast_file);
    ll_utm=ll2utm(ll(:,1:2),origin);
elseif ~exist('ll') && plot_figures==0
    ll=[];
    ll_utm=[];
end

%% Gutenburg-Richter Plots

% Set Histogram Parameters (0.1 Mw Bins)
bins=[floor(min(eqs(:,4))*10)/10:0.1:ceil(max(eqs(:,4))*10)/10];

[n,~,~,btrend,mdl]=gutenberg_richter(eqs(:,4),bins,bvalue,[],MOC,plot_figures);

if isempty(avalue)
    avalue=mdl.Coefficients.Estimate(1);
end

if isempty(bvalue)
    bvalue=-mdl.Coefficients.Estimate(2)
end

%% Calculate displacement of fractures MAY BE WORTH SPLITTING INTO SEPERATE FUNCTION SOON ENOUGH
% https://www.sciencedirect.com/science/article/pii/S0191814119303657
% dmax=0.0337L^1.02
% I am starting to suspect the validity of this law (specifically the use
% of Dmax). maybe try disp=3e-5 L (Rich's paper
% https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/jgrb.50236 or https://www.sciencedirect.com/science/article/pii/S0040195113003776)

Mag=bins; % Get the seismic magnitude
mo=10.^((3/2)*(Mag+6.07)); % Convert to seismic moment
M=mo/(pi*mu*2*scaling_a);
len=nthroot(8*M,scaling_exp+2);
disp=scaling_a.*(len.^scaling_exp); % Use scaling law to work out displacement
mag_window=(2/3)*log10((16/7)*(7/16)*pi*mu*(scaling_a*(vein_lengths.^scaling_exp)/(vein_lengths*0.5))*(0.5*vein_lengths).^3)-6.07;

fprintf('For Fractures %.1f - %.1f m, magnitudes = Mw %.1f to %.1f\n', vein_lengths(1),vein_lengths(2),mag_window(1),mag_window(2))



%% Quadtree subsampling of the data

% Rotation flags
r_flg.norm2strike=1; % Rotate dataset to strike before quadtree
r_flg.norm2dip=0; % Rotate dataset to dip before quadtree
r_flg.strike2norm=0; % Rotate dataset back from strike to normal after quadtree
r_flg.dip2norm=0; % Rotate dataset back from dip to normal after quadtree

[QT,dropped] = quadtree_subsample(eqs(:,1:3),af,bearing,dip,QuadminbinCapacity,r_flg,ll_utm,'shrink',plot_figures,'centeredSquare',50);


%% Fill gaps
fill_gaps=1;
if fill_gaps==1
% figure();QT.plot; hold on; plot(eqs(dropped,1),eqs(dropped,2),'k.'); title('Iteration 1')
[DT,dropped2] = quadtree_subsample(eqs(dropped,1:3),af,bearing,dip,QuadminbinCapacity,r_flg,ll_utm,'shrink',plot_figures,'centeredSquare',50);

QT.Points=[QT.Points;DT.Points];
QT.PointBins=[QT.PointBins;DT.PointBins+QT.BinCount];
QT.BinCount=QT.BinCount+DT.BinCount;
QT.BinBoundaries=[QT.BinBoundaries;DT.BinBoundaries];
QT.BinCorners=cat(3,QT.BinCorners,DT.BinCorners);

eqd=eqs(dropped,1:3);
% figure();QT.plot; hold on; plot(eqd(dropped2,1),eqd(dropped2,2),'k.'); title('Iteration 2')

[DT,dropped3] = quadtree_subsample(eqd(dropped2,1:3),af,bearing,dip,QuadminbinCapacity,r_flg,ll_utm,'shrink',plot_figures,'centeredSquare',50);

QT.Points=[QT.Points;DT.Points];
QT.PointBins=[QT.PointBins;DT.PointBins+QT.BinCount];
QT.BinCount=QT.BinCount+DT.BinCount;
QT.BinBoundaries=[QT.BinBoundaries;DT.BinBoundaries];
QT.BinCorners=cat(3,QT.BinCorners,DT.BinCorners);

% eqd=eqd(dropped2,1:3);
% figure();QT.plot; hold on; plot(eqd(dropped3,1),eqd(dropped3,2),'k.');; title('Iteration 3')
ix=(1:length(dropped2))';
ix=ismember(ix,dropped3);
dropped2(ix==0)=[];
ix=(1:length(dropped))';
ix=ismember(ix,dropped2);
dropped(ix==0)=[];
eqs(dropped,:)=[];
end

%% Determine Upper and Lower Limits of Seismicity

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
% save('depth_file.mat','seis_depths')
% fid=fopen('depth_fileu.txt','w');
% for ii=1:QT.BinCount
% fprintf(fid,'%.0f %.0f %.0f\n',seis_depths(ii,1),seis_depths(ii,2),seis_depths(ii,3));
% end
% fclose(fid);
% fid=fopen('depth_filel.txt','w')
% for ii=1:QT.BinCount
% fprintf(fid,'%.0f %.0f %.0f\n',seis_depths(ii,1),seis_depths(ii,2),seis_depths(ii,4));
% end
% fclose(fid);
%%

exhum=readmatrix('C:\Jacko\NZ 2020\seis2frac\Modeled_exhumation_rates_grid_points.csv');
exhum(:,(1:2))=ll2utm(exhum(:,(1:2)),origin);

[thickness,fracture_time,ex]=exhumation_time(bin_center,exhum,seis_depths,seis_ix,QT,1,af,bearing,origin,ll);

bin_total_fractures=(fracture_time/10).*bin_total;
bin_fracture_density=bin_total_fractures./(bin_area.*thickness);

%%
transect_density=[];

for ii=1:QT.BinCount
    [transect_density(ii)]=calculate_transect(ii,QT,thickness,eqs,fracture_time,len,bins,bvalue,MOC,vein_lengths,vein_thickness,0);
end


%%
[bin_center]=local2llh(bin_center',origin)';
[eqs(:,1:2)]=local2llh(eqs(:,1:2)',origin)';
[af(:,1:2)]=local2llh(af(:,1:2)',origin)';
for ii=1:QT.BinCount
    [QT.BinCorners(:,1:2,ii)]=local2llh(QT.BinCorners(:,1:2,ii)',origin)';
end

%%
upperlim=20;
figure
title('Transect Fracture Density')
hold on
load C:/Jacko/scripts/ScientificColourMaps6/bamako/bamako.mat;colormap(bamako);
% bamako=colormap(jet);
c=colorbar('FontSize',15);
c.Label.String='Fractures Intersected per Meter';
caxis([0 upperlim])


axis equal

for ii=1:QT.BinCount
    try
        fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*(transect_density(ii)/upperlim)),:))
        
    catch
        if ~isnan(fracture_time(ii))
            fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(size(bamako,1),:))
        end
    end
end
try
    plot(af_ll([1 2],1),af_ll([1 2],2),'r');
end

plot(170.14,-43.60,'kp','MarkerFaceColor','g','MarkerEdgeColor','w')
xlim([ 168.1365  172.1054])
ylim([ -44.7979  -42.0618])
plot(ll(:,1),ll(:,2),'k')

%% Figures

if plot_figures==1
    
    % Load Coastlines
    if ~exist('ll')
        [lon,lat]=get_coast(169,173,-45,-42,[],'f',coast_file);
        ll=ll2utm([lon,lat],origin);
    end
    
    %%
    figure;
    figname='Seismicity Distribution';
    set(gcf,'renderer','zbuffer','name',figname); title(figname);
    hold on
    scatter3(eqs(:,1),eqs(:,2),eqs(:,3),10*(eqs(:,4)+1.5),eqs(:,3),'filled'); colormap(flipud(jet)); c=colorbar; c.Label.String='Depth (km)';
    plot(ll(:,1),ll(:,2))
    try
        plot3(af([1 2 4 3 1],1),af([1 2 4 3 1],2),af([1 2 4 3 1],3));
    end
    xlabel('Longitude');ylabel('Latitude');zlabel('Depth');view([-40,10]);%axis equal
    hold off
    
    %%
    f=figure;
    figname='Displacement-Length-Magnitude';
    set(gcf,'renderer','zbuffer','name',figname); title(figname,'FontSize',13);
    hold on
    scatter(log10(len),log10(disp),20,bins,'filled')
    xlabel('Fracture Length (m)','FontSize',13); ylabel('Fracture Displacement (m)','FontSize',13); c=colorbar; c.Label.String='Seismic Magnitude';c.Label.FontSize=13;
    xtickpoints=[-5:1:5]; ytickpoints=[-6:1:4];
    xticks(xtickpoints); xticklabels(10.^xtickpoints);
    yticks(ytickpoints); yticklabels(10.^ytickpoints);
    for ii=1:2:length(xtickpoints); xline(xtickpoints(ii)); end;
    for ii=1:2:length(ytickpoints); yline(ytickpoints(ii)); end;
    xlim([-5 5])
    ylim([-6 2])
    
    %%
    %     figure
    %     scatter(bin_center(:,1),bin_center(:,2),bin_center(:,3),10*bin_volume/mean(bin_volume),log10(bin_total),'filled')
    %     title('Projected Events per grid')
    %     xlabel('Longitude');ylabel('Latitude');zlabel('Depth');%%pbaspect([1 1 1])
    %     hold on
    %     plot(ll(:,1),ll(:,2))
    %     c=colorbar; c.Label.String='Projected Events';
    %     caxis([0 ceil(max(bin_alpha))])
    %     c.Ticks=[0:1:ceil(max(bin_alpha))];
    %     for ii=1:ceil(max(bin_alpha))+1
    %         c.TickLabels{ii}=num2str(10.^c.Ticks(ii));
    %     end
    %
    %     plot3(100,40,3,'rp','MarkerFaceColor','r')
    %%
    figure
    title('Fracture Density in column')
    xlabel('Longitude');ylabel('Latitude');
    hold on
    bamako=colormap(jet);colorbar;
    c=colorbar;
    caxis([0 5e6])
    
    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*(bin_fracture_density(ii)/5e6)),:))
            
        catch
            if ~isnan(fracture_time(ii))
                fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(size(bamako,1),:))
            end
        end
    end
    try
        plot(af([1 2],1),af([1 2],2),'r');
    end
    plot(ll(:,1),ll(:,2),'k')
    %%
    figure
    title('Fracture Total in column')
    xlabel('Longitude');ylabel('Latitude');
    hold on
    bamako=colormap(jet);colorbar; %caxis([0 25])
    c=colorbar;
    caxis([0 1e9])
    
    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*(bin_total_fractures(ii)/1e9)),:))
            
        catch
            if ~isnan(fracture_time(ii))
                fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(size(bamako,1),:))
            end
        end
    end
    try
        plot(af([1 2],1),af([1 2],2),'r');
    end
    plot(ll(:,1),ll(:,2),'k')
end
%%
figure
title('Quad Bin Numbers')
xlabel('Longitude');ylabel('Latitude');
hold on
plot(ll(:,1),ll(:,2))

try
    plot(af_ll([1 2],1),af_ll([1 2],2));
end
axis equal

for ii=1:QT.BinCount
    try
        plot(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii))
        text(bin_center(ii,1),bin_center(ii,2),num2str(ii),'FontSize',5)
    end
end
%%
figure
title('Fracture Transect Density')
xlabel('Longitude');ylabel('Latitude');
hold on
plot(ll(:,1),ll(:,2))

try
    plot(af_ll([1 2],1),af_ll([1 2],2));
end
axis equal

for ii=1:QT.BinCount
    try
        plot(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii))
        text(bin_center(ii,1),bin_center(ii,2),num2str(round(transect_density(ii)*10)/10),'FontSize',5)
    end
end