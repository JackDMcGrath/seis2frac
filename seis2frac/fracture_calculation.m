function [transect_density,volume_change]=fracture_calculation(poly_width,step_distance,step_time,MOC,MOC_ix,bins,D,ray,observation_time,bvalue,len,len_ix,vein_thickness,frac_area,disrate,molar_density)


ray_length=sum(step_distance); % km
fracture_total_time=sum(step_time); % Years
ray_volume=poly_width^2*ray_length; % Volume in which fractures occur
tmp_eq=[];
% volume_change=1;

for ii=1:size(MOC_ix,2)-1 % FOR EACH MAGNITUDE, CALCULATE DENSITY OF EVENTS FOR EACH TIMESTEP
    MOC_eqs(:,ii)=D{ii}(ray); % Density of events at each timestep (f/km^3) 1 cell per ray, 1 col per mag
    MOC_eqs_total=round(mean(MOC_eqs).*ray_volume*fracture_total_time/observation_time); % Multiply mean density by volume and fracture time to get total events > MOC
    tmp_eq=[tmp_eq,linspace(bins(MOC_ix(ii)),bins(MOC_ix(ii)+1),MOC_eqs_total(ii))];
end

if size(find(MOC_eqs_total),2)==0
    transect_density=0;
    fracture_total_volume=0;
    disvol_total=0;
    volume_change=1;
    fprintf(' No events found around raypath\n')
else
    if size(find(MOC_eqs_total),2)==1
        b=bvalue;
        fprintf(' Events only in one magnitude bin. Global b-value used\n')
    else
        b=[];
    end
    [~,~,~,bin_btrend]=gutenberg_richter(tmp_eq',bins,b,[],MOC,0); % Find B-trend of this volume. B-values calculated from the input data
    bin_btrend=bin_btrend(len_ix).^10; % Number of events in each pathpoly for each magnitude, corrected for b-value
    eq_tot=sum(bin_btrend); % Total corrected number of events in entire volume
    fracture_volumes=pi*((0.5*len(len_ix)).^2)*vein_thickness*1e-3; % Volume of each fracture length in m^3
    fracture_total_volume=sum(fracture_volumes.*bin_btrend*1e-6); % Total Fracture volume km^3
    
    %% Work out volume change
    step_time_percent=step_time./fracture_total_time; % Percent of total time spent in each time step
    total_area=sum(bin_btrend.*frac_area); % m^2
    
    % Calculate dissovled volume
    for ii=1:size(ray,1)
    disvol(ii)=dissolved_volume(disrate(ii),step_time_percent(ii)*total_area,step_time(ii),molar_density);
    end
    disvol_total=sum(disvol*1e-9);
    volume_change=1-(disvol_total/ray_volume);
    
    fracture_ratio=fracture_total_volume/(ray_volume*volume_change);
    transect_density=fracture_ratio/(vein_thickness*1e-3);
    
end
fprintf(' Rock Volume: %.1f km^3\n Dissolved Volume: %.3f km^3\n Fracture Volume: %.3f km^3\n Transect Density: %.8f f/m\n Volume Change: %.2f\n', ...
    ray_volume,disvol_total,fracture_total_volume,transect_density,volume_change)
end