mu=32e9; % Shear modulus, Pa
scaling_a=3e-5;
scaling_exp=1;
% Mag=[-1.1:0.1:4.6]; % Get the seismic magnitude
% moment=10.^((3*(Mag+6.07))/2); % Convert to seismic moment
% %     len(ii)=10^(log10(moment/(mu*pi*(scaling_exp+2)))); % Calculate length required TYPO'D
% len=10.^(log10(moment/(mu*pi*scaling_a))/(scaling_exp+2)); % Calculate length required
% disp=scaling_a.*(len.^scaling_exp); % Use scaling law to work out displacement
% 
% stress=((7*pi)/16)*mu*(disp/(len));
% Mo=(16/7).*stress.*(0.5*len).^3;
% Mw=(2/3).*log(Mo)-6.07;
% 
% % plot(Mag,Mw,'.')
% 
L=43.892

M=L^3*(0.5*L)^2

l=10^((log10(M)-2*log10(0.5))/5)
% 
Mag=[-1.1:0.1:4.6]; % Get the seismic magnitude
mo=10.^((3/2)*(Mag+6.07)); % Convert to seismic moment
M=mo/(pi*mu*2*scaling_a);
len=nthroot(8*M,scaling_exp+2);
disp=scaling_a.*(len.^scaling_exp); % Use scaling law to work out displacement

radius=len/2;
stress=((7*pi)/16)*mu*(disp/radius);
Mo=(16/7).*stress.*(radius.^3);
Mw=(2/3).*log10(Mo)-6.07;

plot(Mag,Mw,'.')
hold on
plot(-1:5,-1:5)
Mw-Mag

L=43.892;

M=L^(scaling_exp-1)*(0.5*L)^3;

l=10^((log10(M)-3*log10(0.5))/(scaling_exp+2));


L=18;

a=((scaling_a)*(L^scaling_exp))/(0.5*L)*(0.5*L)^3
a2=((2*scaling_a*(L^scaling_exp))/L)*(0.5*L)^3
a==a2