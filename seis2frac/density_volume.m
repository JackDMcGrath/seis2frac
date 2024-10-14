function [dmap,D]=density_volume(eqs,bins,MOC,MOC_ix)
%% DENSITY_VOLUME
%
% Use the earthquakes to generate the earthquake density volumes

dmapX=[min(round(eqs(:,1))):5:max(round(eqs(:,1)))];
dmapY=[min(round(eqs(:,2))):5:0];
dmapZ=[min(round(eqs(:,3))):1:max(round(eqs(:,3)))];
[x,y,z]=meshgrid(dmapX,dmapY,dmapZ);

for ii=1:(size(MOC_ix,2)-1) % FOR EACH MAGNITUDE BIN, CALCULATE THE INTERPOLATION FUNCTION
    eqs_ix=find(((eqs(:,4)<=bins(MOC_ix(ii)+1))+(eqs(:,4)>=bins(MOC_ix(ii))))==2); % Only look at earthquakes within the magnitude window
    [dmap{ii}]=dataDensity3(eqs(eqs_ix,1),eqs(eqs_ix,2),eqs(eqs_ix,3),size(dmapX,2),size(dmapY,2),size(dmapZ,2),10,[],1); % Calculate density for each magnitude window
    D{ii}=scatteredInterpolant(x(:),y(:),z(:),dmap{ii}(:));
end

[dmap{size(MOC_ix,2)}]=dataDensity3(eqs(eqs(:,4)>MOC,1),eqs(eqs(:,4)>MOC,2),eqs(eqs(:,4)>MOC,3),size(dmapX,2),size(dmapY,2),size(dmapZ,2),10,[],1); % Calculate density for each magnitude window
