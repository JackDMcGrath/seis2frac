%% SEIS2FRAC
% Script for analysing seismicity and relating to fracture propeties
fprintf('This should be branch QUADTREE\n')
%% Setup
% Admin
plot_figures=1;
eq_file='seismic_data.csv';
coast_file='C:\Jacko\scripts\coastlines\gshhg-bin-2.3.7';

% Geometry
gridx=10; % set size of cube or minSize of OcTree
OctbinCapacity=20; % Maximum number of events > MOC in OcTree
QuadminbinCapacity=10;
style='Normal'; % OcTree division method ('Normal' or 'Weighted')
grdShape='Cube'; % Shape of OcTree grids ('Normal' or 'Cube')
dip=60; % set fault dip
fault_base=35; % set depth of fault
search_lims=[250,40,30]; % Area to search in km (x,y,z)

%Rock properties
mu=32e9; % Shear modulus, Pa

% G-R Constants from Michalios 2019
avalue=[]; % If empty, will calculate from dataset (4.28)
bvalue=0.85; % If empty, will calculate from dataset. Value from 9111 events =0.85
MOC=1.1;

% Fault Length-Displacement Scaling Parameters
scaling_a=0.0337;
scaling_exp=1.02;


%% Prepare Data
% Seismic Catalogue
eqs=readmatrix(eq_file); %Import seismicity data
eqs=[eqs(:,8),eqs(:,7),-eqs(:,9),eqs(:,10)]; % Lon, Lat, Depth (-ve), Mag
eqs(:,[5:6])=nan; % Lon, Lat, Depth, Mag, Dist to Fault, Moment

% Fault Limits
af=[169.0844,-43.9019,0;171.9880,-42.4733,0]; % Define ends of fault
if ~isempty(af)
    origin=af(1,1:2); % Set UTM origin as one end of fault
    % Calculation: co-ords of base of fault
    af(:,1:2)=ll2utm(af(:,1:2),origin); % Set to UTM
    bearing=atand((af(2,1)-af(1,1))/(af(2,2)-af(1,2)));
    af(3:4,:)=[af(1,1)+(cosd(bearing)*(fault_base/tand(dip))),af(1,2)-(sind(bearing)*(fault_base/tand(dip))),-fault_base;...
        af(2,1)+(cosd(bearing)*(fault_base/tand(dip))),af(2,2)-(sind(bearing)*(fault_base/tand(dip))),-fault_base];
    % Calculation: Horizontal distance to Fault trace (Map view)
    eqs(:,5)=point_to_line(eqs(:,(1:2)),af(1,:),af(2,:));
else
    origin=mean(eqs(:,[1 2]));
end
% Convert to UTM
eqs(:,(1:2))=ll2utm(eqs(:,(1:2)),origin); % Set to UTM




% Calculation: Seismic Moment of each event
eqs(:,6)=10.^((3/2).*(eqs(:,4)+6.07));

%% Gutenburg-Richter Plots

% Set Histogram Parameters
bins=[floor(min(eqs(:,4))*10)/10:0.1:ceil(max(eqs(:,4))*10)/10];

[n,~,~,btrend,mdl]=gutenberg_richter(eqs(:,4),bins,bvalue,[],MOC,plot_figures);

if isempty(avalue)
    avalue=mdl.Coefficients.Estimate(1);
end

if isempty(bvalue)
    bvalue=-mdl.Coefficients.Estimate(2)
end

missing=10^btrend(1)-size(eqs,1); % Number of events missing from the dataset (ie below MOC)

fprintf('Total number of events in catalogue: %.0f\n',length(eqs))
fprintf('Total number of events "missing": %.0f\n',sum(missing))
fprintf('Total number of events: %.0f\n',length(eqs)+sum(missing))

%% Calculate displacement of fractures
% https://www.sciencedirect.com/science/article/pii/S0191814119303657
% dmax=0.0337L^1.02

for ii=1:length(bins)
    if n(ii)~=0
        Mag=bins(ii);
        moment=10.^((3*(Mag+6.07))/2);
        len(ii)=10^(log10(moment/(mu*pi*(scaling_exp+2))));
        disp(ii)=scaling_a*(len(ii)^scaling_exp);
    end
end


%% OcTree subsampling of the data

[OT,eqs]=octree_subsample(eqs,OctbinCapacity,gridx,plot_figures,MOC,af,style,grdShape);

fprintf('Search Complete: \n %.0f events located into %.0f bins\n',size(eqs,1),OT.BinCount)

%% Bin GR

% bin_data=nan(OT.BinCount,6); % events_measured, total_events, bin_volume, fracture_density ,a-value, events > MOC
bin_events=nan(OT.BinCount,1);
bin_MOC=nan(OT.BinCount,1);
bin_total=nan(OT.BinCount,1);
bin_volume=nan(OT.BinCount,1);
bin_density=nan(OT.BinCount,1);
bin_alpha=nan(OT.BinCount,1);
bin_center=nan(OT.BinCount,3);

for ii=1:OT.BinCount
    pflg=0;
    %     if ii==2 || ii==355 ||ii==360
    %         pflg=1;
    %     end
    
    bin_eqs=eqs(OT.PointBins==ii,:); % Get events in one bin
    
    [bin_n,~,bin_alpha(ii),bin_btrend]=gutenberg_richter(bin_eqs(:,4),bins,bvalue,[],MOC,pflg); % Calculate GR trends for bin
    
    bin_events(ii)=size(bin_eqs,1);
    bin_MOC(ii)=length(find(bin_eqs(:,4)>=MOC));
    bin_total(ii)=10^bin_btrend(1);
    
    bin_volume(ii)=prod([OT.BinBoundaries(ii,4)-OT.BinBoundaries(ii,1), ... % Calculate bin volume (km^3)
        OT.BinBoundaries(ii,5)-OT.BinBoundaries(ii,2), ...
        OT.BinBoundaries(ii,6)-OT.BinBoundaries(ii,3)]);
    bin_density(ii)=(10^bin_btrend(1))/bin_volume(ii); % Fracture density in bin inc. "missing" fractures (fractures/km^3)
    for ii=1:OT.BinCount
        bin_center(ii,:)=mean(OT.BinCorners(:,:,ii));
    end
end


%% Determine Upper and Lower Limits of Seismicity

[QT] = quadtree_subsample(eqs(:,1:3),af(1:2,:),bearing,QuadminbinCapacity);

seis_depths=nan(QT.BinCount,4);

for ii=1:QT.BinCount
    if sum(QT.PointBins==ii)>=QuadminbinCapacity
        seis_depths(ii,1:2)=mean(QT.BinCorners(1:4,:,ii)); % Find bin center
        [seis_depths(ii,3),seis_depths(ii,4)] = upper_lower_seis(eqs(find(QT.PointBins==ii),3));
    end
    quad_center(ii,:)=mean(QT.BinCorners(1:4,1:2,ii));
end

seis_ix=find(~isnan(seis_depths(:,1))); % identify bins that have upper-lower depths
%%

exhum=readmatrix('C:\Jacko\NZ 2020\seis2frac\Modeled_exhumation_rates_grid_points.csv');
exhum(:,(1:2))=ll2utm(exhum(:,(1:2)),origin);

% [thickness,fracture_time,ex]=exhumation_time(bin_center,quad_center,exhum,seis_depths,seis_ix,OT,QT,1,af,origin);

% We now have the seismogenic thickness, exhumation rate and the total
% amount of time availiable for fracturing for the OctTree
plot(thickness,ex,'.')

% Use fracture time and the number of events per cube in 10 years to work
% out total fractures

bin_total_fractures=(fracture_time/10).*bin_total;
bin_fracture_density=bin_total_fractures./bin_volume;

%% Figures

if plot_figures==1
    
    % Load Coastlines
    if ~exist('ll')
        [lon,lat]=get_coast(169,173,-45,-42,[],'f',coast_file);
        ll=ll2utm([lon,lat],origin);
    end
    
    figure;
    figname='Seismicity Distribution';
    set(gcf,'renderer','zbuffer','name',figname); title(figname);
    hold on
    scatter3(eqs(:,1),eqs(:,2),eqs(:,3),5*(eqs(:,4)+1.5),eqs(:,3)); colormap(flipud(jet)); c=colorbar; c.Label.String='Depth (km)';
    plot(ll(:,1),ll(:,2))
    try
        plot3(af([1 2 4 3 1],1),af([1 2 4 3 1],2),af([1 2 4 3 1],3));
    end
    xlabel('Lon');ylabel('Lat');zlabel('Depth'); axis equal;view([-40,10])
    hold off
    
    figure;
    figname='Displacement-Length-Magnitude';
    set(gcf,'renderer','zbuffer','name',figname); title(figname);
    hold on
    scatter(log10(len(n~=0)),log10(disp(n~=0)),20,bins(n~=0),'filled')
    xlabel('Fracture Length (m)'); ylabel('Displacement (m)'); c=colorbar; c.Label.String='Magnitude';
    xtickpoints=[-5:2:5]; ytickpoints=[-6:2:4];
    xticks(xtickpoints); xticklabels(10.^xtickpoints);
    yticks(ytickpoints); yticklabels(10.^ytickpoints);
    for ii=1:length(xtickpoints); xline(xtickpoints(ii)); end;
    for ii=1:length(ytickpoints); yline(ytickpoints(ii)); end;
    
    %%
    figure
    scatter3(bin_center(:,1),bin_center(:,2),bin_center(:,3),10*bin_volume/mean(bin_volume),log10(bin_total),'filled')
    title('Projected Events per grid')
    xlabel('Lon');ylabel('Lat');zlabel('Depth');%%pbaspect([1 1 1])
    hold on
    plot(ll(:,1),ll(:,2))
    c=colorbar; c.Label.String='Projected Events';
    caxis([0 ceil(max(bin_alpha))])
    c.Ticks=[0:1:ceil(max(bin_alpha))];
    for ii=1:ceil(max(bin_alpha))+1
        c.TickLabels{ii}=num2str(10.^c.Ticks(ii));
    end
    
    plot3(100,40,3,'rp','MarkerFaceColor','r')
    %%
    figure
    title('Quad Bin Numbers')
    xlabel('Lon');ylabel('Lat');
    hold on
    plot(ll(:,1),ll(:,2))
    
    try
        plot(af([1 2],1),af([1 2],2));
    end
    axis equal
    
    for ii=1:QT.BinCount
        try
            plot(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii))
            text(quad_center(ii,1),quad_center(ii,2),num2str(ii),'FontSize',5)
        end
    end
end

