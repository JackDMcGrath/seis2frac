%% function data_density_transect_estimation(bins,mag_window,eqs,MOC,vein_lengths,len,ray,poly_width,ftime,observation_time,bvalue,volume_change,vein_thickness)
% DATA_DENSITY_TRANSECT_EXTIMATION
% Script to estimate the number of fractures occurring in a transect, using
% a seismicity density volume

%% Calculate the dissolution rate at each temperature step

for ii=1:size(PT.temp,2)
disrate{ii}=dissolution_rate(PT.temp{ii})
end

fprintf('Seeming error in dissolution rate. Fudge applied (divide LogdisRate by 10\n')
fprintf('Also, this equation is for (relatively) low pressures\n')

% Calculate how many moles of silica are in a cm^3 of quartz
quartz_density=2.62; % g/cm^3
silica_molar_weight=60.0843; % g/mol
molar_density=silica_molar_weight/quartz_density; % moles/cm^3

%% Version using dataDensity

% Identify the bins of interest
mag_ix=find(((bins<=mag_window(2))+(bins>=mag_window(1)))==2); % Only look at earthquakes within the magnitude window
mag_ix=sort([mag_ix, mag_ix(1)-1, mag_ix(end)+1]);
MOC_ix=find(((bins<=max(eqs(:,4)))+(bins>=MOC))==2); % Get bins above MOC
MOC_ix=sort([MOC_ix, MOC_ix(end)+1]);
len_ix=find(((len<vein_lengths(2))+(len>vein_lengths(1)))==2); % Find the bins that are of the correct lengths

dmapX=[min(round(eqs(:,1))):5:max(round(eqs(:,1)))];
dmapY=[min(round(eqs(:,2))):5:0];
dmapZ=[min(round(eqs(:,3))):1:max(round(eqs(:,3)))];
[x,y,z]=meshgrid(dmapX,dmapY,dmapZ);

for ii=1:(size(MOC_ix,2)-1) % FOR EACH MAGNITUDE BIN, CALCULATE THE INTERPOLATION FUNCTION
    %     eqs_ix=find(((eqs(:,4)<=bins(mag_ix(ii)+1))+(eqs(:,4)>=bins(mag_ix(ii))))==2); % Only look at earthquakes within the magnitude window
    eqs_ix=find(((eqs(:,4)<=bins(MOC_ix(ii)+1))+(eqs(:,4)>=bins(MOC_ix(ii))))==2); % Only look at earthquakes within the magnitude window
    [dmap{ii}]=dataDensity3(eqs(eqs_ix,1),eqs(eqs_ix,2),eqs(eqs_ix,3),size(dmapX,2),size(dmapY,2),size(dmapZ,2),10,[],1); % Calculate density for each magnitude window
    D{ii}=scatteredInterpolant(x(:),y(:),z(:),dmap{ii}(:));
end
[dmap{size(MOC_ix,2)}]=dataDensity3(eqs(eqs(:,4)>MOC,1),eqs(eqs(:,4)>MOC,2),eqs(eqs(:,4)>MOC,3),size(dmapX,2),size(dmapY,2),size(dmapZ,2),10,[],1); % Calculate density for each magnitude window

%% Calculate Fracture Time then GR

for jj=1:size(ray,2) % FOR EACH RAY
    rock_volume(jj)=poly_width^2*sum(ray{jj}(ftime{jj},5)); % Volume in which the fractures occur
    tmp_eq=[];
    if sum(ftime{jj})>1
        fprintf('\n*** Ray %.0f *** (%.0f, %.0f)\n',jj,ray{jj}(1,1),ray{jj}(1,2))
        for ii=1:(size(MOC_ix,2)-1) % FOR EACH MAGNITUDE, CALCULATE DENSITY OF EVENTS FOR EACH TIMESTEP
            MOC_eqs{jj}(:,ii)=D{ii}(ray{jj}(ftime{jj},1),ray{jj}(ftime{jj},2),ray{jj}(ftime{jj},3)); % Density of events at each timestep (f/km^3) 1 cell per ray, 1 col per mag
            MOC_eqs_total{jj}=round(mean(MOC_eqs{jj}).*rock_volume(jj)*sum(ray{jj}(ftime{jj},6))/observation_time); % Multiply mean density by volume and fracture time to get total events > MOC
            tmp_eq=[tmp_eq,linspace(bins(MOC_ix(ii)),bins(MOC_ix(ii)+1),MOC_eqs_total{jj}(ii))];
        end
        if size(find(MOC_eqs_total{jj}),2)==0
            transect_density(jj)=0;
            fracture_total_volume(jj)=0;
            fprintf(' No events found around raypath\n')
        else
            if size(find(MOC_eqs_total{jj}),2)==1
                b=bvalue;
                fprintf(' Events only in one magnitude bin. Global b-value used\n')
            else
                b=[];
            end
            [~,~,~,bin_btrend{jj}]=gutenberg_richter(tmp_eq',bins,b,[],MOC,0); % Find B-trend of this volume. B-values calculated from the input data
            bin_btrend{jj}=bin_btrend{jj}(len_ix).^10; % Number of events in each pathpoly for each magnitude, corrected for b-value
            eq_tot(jj)=sum(bin_btrend{jj}); % Total corrected number of events
            fracture_volumes=pi*((0.5*len(len_ix)).^2)*vein_thickness*1e-3; % Volume of each fracture length in m^3
            fracture_total_volume(jj)=sum(fracture_volumes.*bin_btrend{jj}*1e-6); % Total Fracture volume km^3
            fracture_ratio(jj)=fracture_total_volume(jj)/(rock_volume(jj)*volume_change);
            transect_density(jj)=fracture_ratio(jj)/(vein_thickness*1e-3);
        end
        fprintf(' Rock Volume: %.1f km^3\n Fracture Volume: %.3f km^3\n Transect Density: %.8f f/m\n', ...
            rock_volume(jj),fracture_total_volume(jj),transect_density(jj))
    end
end







%% Version using dataDensity

% Identify the bins of interest
mag_ix=find(((bins<=mag_window(2))+(bins>=mag_window(1)))==2); % Only look at earthquakes within the magnitude window
mag_ix=sort([mag_ix, mag_ix(1)-1, mag_ix(end)+1]);
MOC_ix=find(((bins<=max(eqs(:,4)))+(bins>=MOC))==2); % Get bins above MOC
MOC_ix=sort([MOC_ix, MOC_ix(end)+1]);
len_ix=find(((len<vein_lengths(2))+(len>vein_lengths(1)))==2); % Find the bins that are of the correct lengths

dmapX=[min(round(eqs(:,1))):5:max(round(eqs(:,1)))];
dmapY=[min(round(eqs(:,2))):5:0];
dmapZ=[min(round(eqs(:,3))):1:max(round(eqs(:,3)))];
[x,y,z]=meshgrid(dmapX,dmapY,dmapZ);

for ii=1:(size(MOC_ix,2)-1) % FOR EACH MAGNITUDE BIN, CALCULATE THE INTERPOLATION FUNCTION
    %     eqs_ix=find(((eqs(:,4)<=bins(mag_ix(ii)+1))+(eqs(:,4)>=bins(mag_ix(ii))))==2); % Only look at earthquakes within the magnitude window
    eqs_ix=find(((eqs(:,4)<=bins(MOC_ix(ii)+1))+(eqs(:,4)>=bins(MOC_ix(ii))))==2); % Only look at earthquakes within the magnitude window
    [dmap{ii}]=dataDensity3(eqs(eqs_ix,1),eqs(eqs_ix,2),eqs(eqs_ix,3),size(dmapX,2),size(dmapY,2),size(dmapZ,2),10,[],1); % Calculate density for each magnitude window
    D{ii}=scatteredInterpolant(x(:),y(:),z(:),dmap{ii}(:));
end
[dmap{size(MOC_ix,2)}]=dataDensity3(eqs(eqs(:,4)>MOC,1),eqs(eqs(:,4)>MOC,2),eqs(eqs(:,4)>MOC,3),size(dmapX,2),size(dmapY,2),size(dmapZ,2),10,[],1); % Calculate density for each magnitude window

%% Calculate Fracture Time then GR

for jj=1:size(ray,2) % FOR EACH RAY
    rock_volume(jj)=poly_width^2*sum(ray{jj}(ftime{jj},5)); % Volume in which the fractures occur
    tmp_eq=[];
    if sum(ftime{jj})>1
        fprintf('\n*** Ray %.0f *** (%.0f, %.0f)\n',jj,ray{jj}(1,1),ray{jj}(1,2))
        for ii=1:(size(MOC_ix,2)-1) % FOR EACH MAGNITUDE, CALCULATE DENSITY OF EVENTS FOR EACH TIMESTEP
            MOC_eqs{jj}(:,ii)=D{ii}(ray{jj}(ftime{jj},1),ray{jj}(ftime{jj},2),ray{jj}(ftime{jj},3)); % Density of events at each timestep (f/km^3) 1 cell per ray, 1 col per mag
            MOC_eqs_total{jj}=round(mean(MOC_eqs{jj}).*rock_volume(jj)*sum(ray{jj}(ftime{jj},6))/observation_time); % Multiply mean density by volume and fracture time to get total events > MOC
            tmp_eq=[tmp_eq,linspace(bins(MOC_ix(ii)),bins(MOC_ix(ii)+1),MOC_eqs_total{jj}(ii))];
        end
        if size(find(MOC_eqs_total{jj}),2)==0
            transect_density(jj)=0;
            fracture_total_volume(jj)=0;
            fprintf(' No events found around raypath\n')
        else
            if size(find(MOC_eqs_total{jj}),2)==1
                b=bvalue;
                fprintf(' Events only in one magnitude bin. Global b-value used\n')
            else
                b=[];
            end
            [~,~,~,bin_btrend{jj}]=gutenberg_richter(tmp_eq',bins,b,[],MOC,0); % Find B-trend of this volume. B-values calculated from the input data
            bin_btrend{jj}=bin_btrend{jj}(len_ix).^10; % Number of events in each pathpoly for each magnitude, corrected for b-value
            eq_tot(jj)=sum(bin_btrend{jj}); % Total corrected number of events
            fracture_volumes=pi*((0.5*len(len_ix)).^2)*vein_thickness*1e-3; % Volume of each fracture length in m^3
            fracture_total_volume(jj)=sum(fracture_volumes.*bin_btrend{jj}*1e-6); % Total Fracture volume km^3
            fracture_ratio(jj)=fracture_total_volume(jj)/(rock_volume(jj)*volume_change);
            transect_density(jj)=fracture_ratio(jj)/(vein_thickness*1e-3);
        end
        fprintf(' Rock Volume: %.1f km^3\n Fracture Volume: %.3f km^3\n Transect Density: %.8f f/m\n', ...
            rock_volume(jj),fracture_total_volume(jj),transect_density(jj))
    end
end