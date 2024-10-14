function [vx,vy,vz]=velocity_fields(bearing,origin,xx,yy,gpsfile,exhum)
%origin=[168.685,-44.115];

gps=readmatrix(gpsfile);
[in]=inpolygon(gps(:,1),gps(:,2),[169,169,173.5,173.5,169],[-45,-42,-42,-45,-45]);
fprintf('To do: Unhardcode gps crop [velocity_fields.m]\n')
fprintf('To do: Work out a way to have velocity at the fault equal to plate rate\n')
fprintf('To do: Should the x-vel be 0 (exhumation location moves with rocks)?\n')
gps=gps(in,:);
gps(:,3:4)=[gps(:,3)*cosd(90-bearing)+gps(:,4)*sind(90-bearing),-gps(:,3)*sind(90-bearing)-gps(:,4)*cosd(90-bearing)];
gps(:,1:2)=ll2utm(gps(:,(1:2)),origin(:,1:2));
[gpsrot]=point_rotation(gps(:,1:2),90-bearing,'strike');
gps(:,1:2)=gpsrot(:,1:2);
gps_xinterp=scatteredInterpolant(gps(:,1),gps(:,2),gps(:,3));
gps_vx=gps_xinterp(xx(:),yy(:))*1e-6*0; % Set to 0 if implying that the zone of exhumation is moving with the rocks. Set to 1 if the exhumation location is fixed and the rocks move through it
gps_yinterp=scatteredInterpolant(gps(:,1),gps(:,2),gps(:,4));
gps_vy=gps_yinterp(xx(:),yy(:))*1e-6;

[exrot]=point_rotation(exhum(:,1:2),90-bearing,'strike');
exhum(:,1:2)=exrot(:,1:2);
exhum_interp=scatteredInterpolant(exhum(:,1),exhum(:,2),exhum(:,3));
exhum_interp.ExtrapolationMethod='nearest'; % options nearest or linear
ex_vz=exhum_interp(xx(:),yy(:))*1e-6;

% dip_slip=readmatrix('C:\Jacko\NZ 2020\seis2frac3\dip_slip_rates.csv');
% dip_slip(:,(2:3))=ll2utm(dip_slip(:,(2:3)),origin(:,1:2));
% [diprot]=point_rotation(dip_slip(:,2:3),90-bearing,'strike');
% dip_slip(:,2:3)=diprot(:,1:2);
% dip_slip(:,3)=0; % Put all slip rates on the "fault"
% sliprate=dip_slip(~isnan(dip_slip(:,4)),[2 3 4]); % Slip rate along the fault
% dipv=dip_slip(~isnan(dip_slip(:,5)),[2 3 5]);
% dipn=dipv;
% dipv(:,3)=sind(60)*dipv(:,3); % Uplift due to fault slip
% dipn(:,3)=cos(60)*dipn(:,3); % Convergence due to fault slip
% slip_distance=60; % Distance at which fault velocity contribution is removed
% 
% sliprate_grid=zeros(size(xx));
% dipv_grid=zeros(size(xx));
% dipn_grid=zeros(size(xx));
% 
% f_ix=find(yy(:,1)==0); % Find row equating to fault
% sliprate_grid(f_ix,:)=interp1(sliprate(:,1),sliprate(:,3),xx(f_ix,:)); % Linear interp to grid
% [mdl]=fitlm(xx(1,:),sliprate_grid(f_ix,:)); % Linear model though these (given overall trend is vaguely linear
% sliprate_grid(f_ix,:)=mdl.Coefficients.Estimate(1)+xx(1,:)*mdl.Coefficients.Estimate(2);
% dipv_grid(f_ix,:)=interp1(dipv(:,1),dipv(:,3),xx(f_ix,:),'nearest');
% dipn_grid(f_ix,:)=interp1(dipn(:,1),dipn(:,3),xx(f_ix,:),'nearest');
% 
% for ii=1:size(xx,2)
%     sliprate_grid(:,ii)=(2*sliprate_grid(f_ix,ii)/pi)*atan(abs(yy(:,ii))*1e3/15e3);
%     dipv_grid(:,ii)=(2*dipv_grid(f_ix,ii)/pi)*atan(yy(:,ii)*1e3/15e3);
%     dipn_grid(:,ii)=(2*dipn_grid(f_ix,ii)/pi)*atan(yy(:,ii)*1e3/15e3);
% end
% 
% sliprate_grid=sliprate_grid-mdl.Coefficients.Estimate(1)+xx(1,:)*mdl.Coefficients.Estimate(2);
% dipv_grid=dipv_grid+interp1(dipv(:,1),dipv(:,3),xx(f_ix,:),'nearest');
% dipn_grid=dipn_grid+interp1(dipn(:,1),dipn(:,3),xx(f_ix,:),'nearest');


clear gpsrot exrot diprot

vx=reshape(gps_vx,size(xx));%*0+-35e-6;
vy=reshape(gps_vy,size(xx));%*0+10e-6;
vz=reshape(ex_vz,size(xx));%*0;

end