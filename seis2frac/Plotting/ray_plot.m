function ray_plot(ray,xy,eqs,MOC,fault,coast,dmap,dmapX,dmapY,dmapZ)

raylim=[];
for ii=1:size(ray,2)
    raylim=[raylim;ray{ii}(:,1:3)];
end

for ii=1:size(ray,2)
rayplot{ii}=ray{ii}(:,1:3);
end
figure
streamline(rayplot)
hold on
plot(xy(:,1),xy(:,2),'r.')
plot3(eqs(eqs(:,4)>=MOC,1),eqs(eqs(:,4)>=MOC,2),eqs(eqs(:,4)>=MOC,3),'r.')
fill3(fault([1 2 4 3 1],1),fault([1 2 4 3 1],2),fault([1 2 4 3 1],3),'g')
plot(coast(:,1),coast(:,2),'k')
axis equal
xlim([min(raylim(:,1))-25 max(raylim(:,1))+25])
ylim([min(raylim(:,2))-10 max(raylim(:,2))+20])
% view([-45,15])
view([-90 1])
ylabel('Cross Strike Distance (km)')
zlabel('Depth (km)')

if exist('oversize','var') && ~isempty(oversize)
    slice(dmapX,dmapY,dmapZ,dmap{size(dmap,2)},[xy(oversize,1)],[],[])
    colorbar
end
