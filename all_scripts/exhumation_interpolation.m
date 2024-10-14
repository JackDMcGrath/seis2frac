exhum=readmatrix('C:\Jacko\NZ 2020\seis2frac3\Modeled_exhumation_rates_grid_points.csv');
[exhum(:,1:2)]=ll2utm(exhum(:,1:2),mean(exhum(:,1:2)));
[exrot]=point_rotation(exhum(:,1:2),90-55,'strike');
exhum=[exrot(:,1:2),exhum(:,3)];
ex=exhum(12:end,:);
% ex=exhum;
n=20;

for ii=1:n
    try
        p=polyfitn(ex(:,1:2),ex(:,3),ii);
        zex=polyvaln(p,ex(:,1:2));
        error(ii)=rms(ex(:,3)-zex);
    catch
        error(ii)=nan;
    end
end


n=find(error==min(error));
p=polyfitn(ex(:,1:2),ex(:,3),n);
[x,y]=meshgrid(-100:2:100,-35:1:35);
z=polyvaln(p,[x(:),y(:)]);

scale=floor(log(abs(range(z)))./log(10))-1;
offset=range(z);



surf(x,y,reshape(z,size(x)));
hold on
plot3(ex(:,1),ex(:,2),ex(:,3),'o','MarkerFaceColor','b')
plot3(exhum(:,1),exhum(:,2),exhum(:,3)*(10^scale)*5-offset,'o','MarkerFaceColor','y')
hold off
xlabel('X');ylabel('Y');zlabel('Exhumation Rate (mm/yr)');
legend('Interpolated Surface\newline','Actual Data Points\newline',['Scaled Data Points\newline(for what it should look like)'])
%%
z=polyvaln(p,[0,-160]);
t=datetime('now');
t=cellstr(t+minutes((365.25*24*60)/(z*1e-6/384400)));

fprintf('Christchurch is going up at a rate of %.1e mm/yr, or %.f million km/yr, or %.f ld/yr\nChristchurch will have landed on the Moon in %.1f minutes at %s\n',...
    z,z*1e-12,z*1e-6/384400,(365.25*24*60)/(z*1e-6/384400),t{1}(end-7:end))
