eq_file='Data/seismic_data.csv';
eqs=readmatrix(eq_file); %Import seismicity data
eqs=[eqs(:,8),eqs(:,7),eqs(:,9)]; % Lon, Lat, Depth (-ve), Mag
origin=[168.685,-44.115];
eqs(:,(1:2))=ll2utm(eqs(:,(1:2)),origin); % Set to UTM
[eqs]=point_rotation(eqs,90-60,'strike');

[dmap]=dataDensity(eqs(:,1),eqs(:,2),301,140,[0 300 -80 60]);

figure();
imagesc([0 300],[-80 60],dmap)
hold on
axis xy
plot(eqs(:,1),eqs(:,2),'k.')


xrange=[0:2:300];
yrange=[-80:2:60];
zrange=[0:2:30];
ddlims=[xrange(1) xrange(end) yrange(1) yrange(end) zrange(1) zrange(end)];

[x,y,z]=meshgrid(xrange,yrange,zrange);

[dmap]=dataDensity3(eqs(:,1),eqs(:,2),eqs(:,3),size(x,2),size(x,1),size(x,3),100,ddlims);

figure()
s=slice(x,y,z,dmap,[0 100 200 300],[0 -25],[10]);
set(s,'edgecolor','none');
colorbar
load C:/Jacko/scripts/ScientificColourMaps6/batlow/batlow.mat
colormap(batlow)
