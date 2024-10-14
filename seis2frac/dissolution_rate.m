function [disRate,LogdisRate]=dissolution_rate(temperature)
%% DISSOLUTION RATE
%
% Calculate dissolution rate at various temperatures using relationship of
% Zhang et al (2015 https://www.sciencedirect.com/science/article/pii/S0896844615000716?via%3Dihub
%
% Input: Temperature (degrees celcius)
% Output: Dissolution Rate (mols/m^2/s^1)
% To do: Include variable pressure
%      : Include method to vary constants

temperature=temperature+272.15; % Convert to Kelvin
LogdisRate=(-1.6478-(13573.1./temperature))/10 -2; % I beleive this to be wrong. Divide by 10 and substract 2 is a fudge
disRate=10.^LogdisRate; % mol/m^2/s^1

% figure
% plot(1./temperature,(-1.6478-(13573.1./temperature))/10-2,'k') %% 23MPa (25-374
% hold on
% plot(1./temperature,(-1.325-(14124.0./temperature))/10-2,'r--') %% 33MPa (50-374)
% data=readmatrix('data2.csv');
% plot(1./(data(find(data(:,2)==23),1)+272.15),data(find(data(:,2)==23),3),'sk','MarkerFaceColor','k') %% 33MPa
% plot(1./(data(find(data(:,2)==33),1)+272.15),data(find(data(:,2)==33),3),'or','MarkerFaceColor','r') %% 33MPa
% xlabel('1/T(K)')
% ylabel('Log k+ (mol/m^2/s^1)')
% ylim([-9.25 -2])
% xlim([0.00125 0.00375])
