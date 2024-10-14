function [transect_density]=calculate_transect(bin,QT,thickness,eqs,fracture_time,len,bins,bvalue,MOC,vein_lengths,vein_thickness,pflg)

%% use len(28:34) for only fractures between 1 and 10 m

width=sqrt((QT.BinCorners(1,1,bin)-QT.BinCorners(2,1,bin))^2+(QT.BinCorners(1,2,bin)-QT.BinCorners(2,2,bin))^2);
hgt=thickness(bin);
vol=width*width*hgt;


width_range=0:0.001:width;
width_range=width_range(randperm(length(width_range)));
hgt_range=0:0.001:hgt;
hgt_range=hgt_range(randperm(length(hgt_range)));


bin_eqs=eqs(QT.PointBins==bin,:);
[~,~,~,bin_btrend]=gutenberg_richter(bin_eqs(:,4),bins,bvalue,[],MOC,pflg);

bin_btrend=(10.^bin_btrend)*(fracture_time(bin)/10); % Number of events of each magnitude in the bin across all exhumation time

min_len_ix=min(find(len>=vein_lengths(1)));
max_len_ix=max(find(len<=vein_lengths(2)));
len=len(min_len_ix:max_len_ix);
bin_btrend=bin_btrend(min_len_ix:max_len_ix);

%% Assume every fracture is 1mm wide

fracture_volumes=pi*(len.^2)*vein_thickness*1e-3; % m^2 - set thickness to 1mm

fracture_total_volume=fracture_volumes.*bin_btrend*1e-6; % km^2

fracture_ratio=sum(fracture_total_volume)/vol;

transect_density=1*fracture_ratio/(vein_thickness*1e-3);
end