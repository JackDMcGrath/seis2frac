%% function data_density_transect_estimation(bins,mag_window,eqs,MOC,vein_lengths,len,ray,poly_width,ftime,observation_time,bvalue,volume_change,vein_thickness)
% DATA_DENSITY_TRANSECT_EXTIMATION
% Script to estimate the number of fractures occurring in a transect, using
% a seismicity density volume

%% Calculate the dissolution rate at each temperature step

for ii=1:size(PT.temp,2)
disrate{ii}=dissolution_rate(PT.temp{ii});
end

fprintf('Seeming error in dissolution rate. Fudge applied (divide LogdisRate by 10\n')
fprintf('Also, this equation is for (relatively) low pressures\n')

% Calculate how many moles of silica are in a cm^3 of quartz
quartz_density=2.62; % g/cm^3
silica_molar_weight=60.0843; % g/mol
molar_density=silica_molar_weight/quartz_density; % moles/cm^3

%% Create Data Desnity Volumes

[MOC_ix,len_ix]=identify_fracture_sizes(eqs,bins,mag_window,MOC,len,vein_lengths);
frac_surface_area=pi.*((len(len_ix)*0.5).^2).*2; % Surface area of each fracture length (*2 for left and right hand side) 

[dmap,D]=density_volume(eqs,bins,MOC,MOC_ix);

%% Calculate Fracture Time then GR

transect_density=zeros(size(ray,2),1);

for jj=1:size(ray,2) % FOR EACH RAY
     step_distance=ray{jj}(ftime{jj},5); % Distance while fracturing
     step_time=ray{jj}(ftime{jj},6); % Time while fracturing
     
         if sum(ftime{jj})>1 % If ray enters fracture zone
       
        fprintf('\n*** Ray %.0f *** (%.0f, %.0f)\n',jj,ray{jj}(1,1),ray{jj}(1,2)) % Print Ray ID
        
     [transect_density(jj),disvol_total(jj)]=fracture_calculation(poly_width,step_distance,step_time,MOC,MOC_ix,bins,D,ray{jj}(ftime{jj},1:3),...
         observation_time,bvalue,len,len_ix,vein_thickness,frac_surface_area,disrate{jj}(ftime{jj}),molar_density);
     
    end
end





