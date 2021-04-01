%% SEIS2FRAC
% Script for analysing seismicity and relating to fracture propeties

%% Setup
% Admin
plot_figures=0;
eq_file='seismic_data.csv';
coast_file='C:\Jacko\scripts\coastlines\gshhg-bin-2.3.7';

% Geometry
gridx=1; % set size of cube or minSize of OcTree
binCapacity=20; % Maximum number of events > MOC in OcTree
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
eqs=[eqs(:,8),eqs(:,7),-eqs(:,9),eqs(:,10)]; % Lon, Lat, Depth, Mag
eqs(:,[5:6])=nan; % Lon, Lat, Depth, Mag, Dist to Fault, Moment

% Fault Limits
af=[169.0844,-43.9019;171.9880,-42.4733]; % Define ends of fault
origin=af(1,:); % Set UTM origin as one end of fault

% Convert to UTM
af=ll2utm([169.0844,-43.9019;171.9880,-42.4733],origin); % Set to UTM
af(:,3)=0;
eqs(:,(1:2))=ll2utm(eqs(:,(1:2)),origin); % Set to UTM


% Calculation: co-ords of base of fault
bearing=atand((af(2,1)-af(1,1))/(af(2,2)-af(1,2)));
af(3:4,:)=[af(1,1)+(cosd(bearing)*(fault_base/tand(bearing))),af(1,2)-(sind(bearing)*(fault_base/tand(bearing))),35;...
    af(2,1)+(cosd(bearing)*(fault_base/tand(bearing))),af(2,2)-(sind(bearing)*(fault_base/tand(bearing))),35];

% Calculation: Horizontal distance to Fault trace (Map view)
eqs(:,5)=point_to_line(eqs(:,(1:2)),af(1,:),af(2,:));

% Calculation: Seismic Moment of each event
eqs(:,6)=10.^((3/2).*(eqs(:,4)+6.07));

%% Gutenburg-Richter Plots

% Set Histogram Parameters
bins=[floor(min(eqs(:,4))*10)/10:0.1:ceil(max(eqs(:,4))*10)/10];

[n,~,~,btrend,mdl]=gutenberg_richter(eqs(:,4),bins,bvalue,[],MOC,plot_figures);

if isempty(avalue)
    avalue=mdl.Coefficients.Estimate(1)
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

[OT,eqs]=octree_subsample(eqs,binCapacity,gridx,plot_figures,MOC,af,style,grdShape);

fprintf('Search Complete: \n %.0f events located into %.0f bins\n',size(eqs,1),OT.BinCount)

%% Bin GR
% Work on bin 451 (bin with most (29) events > MOC when grix=1, binCap=5)

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
    if ii==1000
        pflg=1;
    end
    
bin_eqs=eqs(OT.PointBins==ii,:); % Get events in one bin

[bin_n,~,bin_alpha(ii),bin_btrend]=gutenberg_richter(bin_eqs(:,4),bins,bvalue,[],MOC,pflg); % Calculate GR trends for bin

bin_events(ii)=size(bin_eqs,1);
bin_MOC(ii)=length(find(bin_eqs(:,4)>=MOC));
bin_total(ii)=10^bin_btrend(1);

bin_volume(ii)=prod([OT.BinBoundaries(ii,4)-OT.BinBoundaries(ii,1), ... % Calculate bin volume (km^3)
    OT.BinBoundaries(ii,5)-OT.BinBoundaries(ii,2), ...
    OT.BinBoundaries(ii,6)-OT.BinBoundaries(ii,3)]);
bin_density(ii)=(10^bin_btrend(1))/bin_volume(ii); % Fracture density in bin inc. "missing" fractures (fractures/km^3)
bin_center(ii,:)= mean([OT.BinBoundaries(ii,[1:3]);OT.BinBoundaries(ii,[4:6])]); % Calculate bin centers

end
%%
 figure
 scatter3(bin_center(:,1),bin_center(:,2),bin_center(:,3),bin_total,bin_density,'filled')
%  title('Cumulative Seismic Moment per grid')
 xlabel('Lon');ylabel('Lat');zlabel('Depth');%%pbaspect([1 1 1])
 hold on
 plot(ll(:,1),ll(:,2))
 c=colorbar; c.Label.String='Density';




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
    plot3(af([1 2 4 3 1],1),af([1 2 4 3 1],2),-af([1 2 4 3 1],3));
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
        
end

