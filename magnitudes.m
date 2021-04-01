%% MAGNITUDES
% Script for calculating magnitudes of fractures based on a range
% displacements and fracture lengths

%% From clarke 66 (convert Mb to Pa)
mu= [0.39]*1e11; % [shale,shist]
p=0.181;
lame=mu*((2*p)/(1-2*p));

fault_type='c'; %c , s, d


disp=[0.001,0.1];
a=[1:0.1:15];

dstress=nan(length(disp),length(a));

Mo=nan(length(disp),length(a));

Mw=nan(length(disp),length(a));

if strcmpi(fault_type,'c')
for ii=1:length(disp) % Circular fracture
dstress(ii,:)=((7*pi)/16).*mu.*(disp(ii)./a);
Mo(ii,:)=(16/7).*dstress(ii,:).*(a.^3);
Mw(ii,:)=(2/3)*log10(Mo(ii,:))-6.07;
end
elseif strcmpi(fault_type,'s')
for ii=1:length(disp) % Strike-Slip fracture
dstress(ii,:)=(2/pi).*mu.*(disp(ii)./a);
Mo(ii,:)=(pi/2).*dstress(ii,:).*(a.^2)*disp(ii);
Mw(ii,:)=(2/3)*log10(Mo(ii,:))-6.07;
end
elseif strcmpi(fault_type,'d')
for ii=1:length(disp) % Dip fracture
dstress(ii,:)= ((4*(lame+mu))./(pi*(lame+2*mu)))*mu*(disp(ii)./a);
Mo(ii,:)=((pi*(lame+2*mu))/(4*(lame+mu))).*dstress(ii,:).*(a.^2)*disp(ii);
Mw(ii,:)=(2/3)*log10(Mo(ii,:))-6.07;
end
end
%%

figure()
hold on
for ii=1:length(disp)
    plot(a,Mw(ii,:),'.')
end
if strcmpi(fault_type,'c')
title(['Circular Fault, Mu = ' num2str(mu*1e-9) ' GPa'])
elseif strcmpi(fault_type,'s')
    title('Strike-Slip')
elseif strcmpi(fault_type,'d')
    title('Dip-Slip Fault')
end
xlabel('Fracture Length (m)')
ylabel('Moment Magnitude')
l=legend('1 mm','1 cm','Location','SouthEast');
title(l,'Displacements')
%%
% figure
% hold on
% xlabel('Fracture Length (m)')
% ylabel('Stress Drop (GPa)')
% 
% plot(a,dstress(1,:)*1e-9,'r.')
% plot(a,dstress(2,:)*1e-9,'b.')
% 
% yyaxis right
% ylabel('Seismic Moment')
% 
% plot(a,Mo(1,:),'rd')
% plot(a,Mo(2,:),'bd')
% 
% legend('Shale Stress Drop','Schist Stress Drop','Shale Moment','Schist Moment')


% figure
% hold on
% xlabel('Fracture Length (m)')
% ylabel('Stress Drop (GPa)')
% 
% plot(a,dstress(1,:)*1e-9,'r.')
% plot(a,dstress(2,:)*1e-9,'b.')
% 
% yyaxis right
% ylabel('Moment Magnitude')
% 
% plot(a,Mw(1,:),'rd')
% plot(a,Mw(2,:),'bd')
% 
% legend('Shale Stress Drop','Schist Stress Drop','Shale Mag','Schist Mag')
% 
% eq=readmatrix('sup1.csv');
% figure
% histogram(eq(:,10));
% xlabel('Magnitude')
% ylabel('Frequency')
