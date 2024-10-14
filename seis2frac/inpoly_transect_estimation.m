%% Version using inpolygon
for ii=1:size(raystart,1)
    shp{ii}=alphaShape(pathpoly{ii}([1:2,5:end],1),pathpoly{ii}([1:2,5:end],2),pathpoly{ii}([1:2,5:end],3)); % Build alphaShape around raypath
    rock_volume(ii)=poly_width^2*sum(ray{ii}(ftime{ii},5)); % Volume in which the fractures occur
    if size(shp{ii}.Points,1)>4 && rock_volume(ii)~=0
        fprintf('\n*** Ray %.0f ***\n',ii)
        try
            [eqsin_all{ii},~]=inShape(shp{ii},eqs(:,1),eqs(:,2),eqs(:,3)); % Find all events within the alphaShape
            eqsin{ii}=logical(eqsin_all{ii}.*eqMOC);
        catch
            eqsin{ii}=logical(zeros(size(eqs,1),1));
            fprintf('AlphaShape %.0f is empty\n',ii)
        end
        
        if sum(eqsin{ii})>1
            [~,~,~,bin_btrend{ii}]=gutenberg_richter(eqs(eqsin{ii},4),bins,[],[],MOC,0); % Find B-trend of this volume. B-values calculated from the input data
            
            len_ix=find(((len<vein_lengths(2))+(len>vein_lengths(1)))==2); % Find the bins that are of the correct lengths
            bin_btrend{ii}=bin_btrend{ii}(len_ix);
            
            eqs_total{ii}=bin_btrend{ii}.^10; % Number of events in each pathpoly for each magnitude, corrected for b-value
            eq_tot(ii)=sum(eqs_total{ii}); % Total corrected number of events
            total_events{ii}=eq_tot(ii)*sum(ray{ii}(ftime{ii},6))/observation_time; % Total number of events that can occur over all exhumation
            fracture_volumes=pi*((0.5*len(len_ix)).^2)*vein_thickness*1e-3; % Volume of each fracture length in m^3
            fracture_total_volume(ii)=sum(fracture_volumes.*eqs_total{ii}*1e-6); % Total Fracture volume km^3
            fracture_ratio(ii)=fracture_total_volume(ii)/(rock_volume(ii)*volume_change);
            
            transect_density(ii)=fracture_ratio(ii)/(vein_thickness*1e-3);
            fprintf(' Rock Volume: %.1f km^3\n Fracture Volume: %.3f km^3\n Observed Events (all): %.0f \n Observed Events (>MOC): %.0f \n B-total Events: %.0f\n Total Events: %.0f\n Transect Density: %.8f f/m\n', ...
                rock_volume(ii),fracture_total_volume(ii),sum(eqsin_all{ii}),sum(eqsin{ii}),eq_tot(ii),total_events{ii},transect_density(ii))
            
        end
        
    end
end

%
load('Data/bamako.mat');
maxd=ceil(max(transect_density(~isinf(transect_density)))*1000)/1000;
figure
hold on
for ii=1:size(ray,2)
    fill(pathpoly{ii}(1:4,1),pathpoly{ii}(1:4,2),bamako(round((transect_density(ii)/maxd)*255)+1,:))
end
colormap(bamako)
colorbar
caxis([0 maxd])
%
for ii=1:size(ray,2)
rayplot{ii}=ray{ii}(:,1:3);
end
figure
streamline(rayplot)
hold on
plot(raystart(:,1),raystart(:,2),'r.')
plot3(eqs(:,1),eqs(:,2),eqs(:,3),'k.')
for ii=1:size(ray,2)
    plot3(eqs(eqsin_all{ii},1),eqs(eqsin_all{ii},2),eqs(eqsin_all{ii},3),'r.')
    plot3(eqs(eqsin{ii},1),eqs(eqsin{ii},2),eqs(eqsin{ii},3),'b.')
%     plot(shp{ii})
end
fill3(af_rot([1 2 4 3 1],1),af_rot([1 2 4 3 1],2),af_rot([1 2 4 3 1],3),'g')