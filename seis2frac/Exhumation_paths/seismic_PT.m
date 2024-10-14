function [P,temp,F]=seismic_PT(seis_depths,pt)

% Function to calculate the pressure and temperature from seismic depths

seis_depths(isnan(seis_depths(:,3)),:)=[]; %Remove anywhere that doesn't have a depth
 
% upper_T=griddata(seis_depths(:,1),seis_depths(:,2),seis_depths(:,3),pt(:,1),pt(:,2));
% lower_T=griddata(seis_depths(:,1),seis_depths(:,2),seis_depths(:,4),pt(:,1),pt(:,2));

upper=scatteredInterpolant(seis_depths(:,1),seis_depths(:,2),seis_depths(:,3)); % Object to calculate upper limit at each point
lower=scatteredInterpolant(seis_depths(:,1),seis_depths(:,2),seis_depths(:,4)); % Object to calculate lower limit at each point

upper=upper(pt(:,1),pt(:,2)); % Find upper lower limits
lower=lower(pt(:,1),pt(:,2));

% ratio=zeros(length(pt(:,1)),1);

ix=find(pt(:,3)>lower); % Find all points above base of seismicity

temp(ix)=(420./lower(ix)).*pt(ix,3); % Taking 420 as temp of seismic base, calculate geotherm then multiply by depth for all seismic points

ix=find(pt(:,3)<lower); % Find Points beneath seismicity
 
temp(ix)=(420+((550-420)./(lower(ix)+35)).*(lower(ix)-pt(ix,3)))'; % Calculate geotherm beneath sesimicty to 35km, and find temp

P=-pt(:,3)*1000*2800*9.81*1e-6; % Simple pressure calculation (MPa)

F=pt(:,3)>lower;

F=F+(pt(:,3)<upper);

F=(F==2); % Fracture timesteps

end