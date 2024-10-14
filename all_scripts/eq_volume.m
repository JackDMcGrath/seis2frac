eq_file='Data/seismic_data.csv';
eqs=readmatrix(eq_file); %Import seismicity data
eqs=[eqs(:,8),eqs(:,7),-eqs(:,9),eqs(:,10)]; % Lon, Lat, Depth (-ve), Mag
origin=[168.685,-44.115,0];
eqs(:,(1:2))=ll2utm(eqs(:,(1:2)),origin(1:2));
eqs(eqs(:,3)<-30,:)=[];

eqs=round(eqs(:,[1:3])); %Round all eqs to nearest km

[uv,uix,uc]=unique(eqs(:,[1:3]),'rows'); % Unique row, unique ix, unique count

X=floor(min(eqs(:,1))):1:ceil(max(eqs(:,1)));
Y=floor(min(eqs(:,2))):1:ceil(max(eqs(:,2)));
Z=floor(min(eqs(:,3))):1:ceil(max(eqs(:,3)));

eq_ref=eqs(:,1:3)-[X(1),Y(1),Z(1)]+1;

vol=zeros([length(Y),length(X),length(Z)]);

for ii=1:length(uix)
    vol(eq_ref(uix(ii),2),eq_ref(uix(ii),1),eq_ref(uix(ii),3))=sum(uc==ii);
end

% vlim=volumebounds(X,Y,Z,vol);
% figure
% s=slice(X,Y,Z,vol,175,[],-10);
% set(s,'edgecolor','none');
% hold on
% plot3(eqs(:,1),eqs(:,2),eqs(:,3),'k.')
% caxis([0 10])

data=smooth3(vol,'box',3);
figure
patch(isocaps(data,.5),...
   'FaceColor','interp','EdgeColor','none');
p1 = patch(isosurface(data,.5),...
   'FaceColor','blue','EdgeColor','none');
isonormals(data,p1);
view(3); 
axis vis3d tight
camlight left
lighting gouraud