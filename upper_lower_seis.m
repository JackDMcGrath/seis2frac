function [lowerseis,upperseis] = upper_lower_seis(depths)
%% UPPER_LOWER_SEIS
% Script to take the earthquake depths and assign an upper and lower limit
% of seismicity based off the 10-90% limits of the cumulative frequencies

plot_flg=0;

neqs=length(depths); % number of events
rng=floor(neqs*0.1):ceil(neqs*0.9); % indexes of 10-90% cumulative frequency
depths=sort(depths); % Sort events into ascending depth order

mdl=fitlm(rng,depths(rng)); % Fit linear model to depths

lowerseis=mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*0;
upperseis=mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*neqs;
%%
if plot_flg==1
figure;
hold on;
plot(1:neqs,depths,'.')
plot(rng,depths(rng),'.')
plot(rng,mdl.Coefficients.Estimate(1)+mdl.Coefficients.Estimate(2)*rng,'g.')
plot(0,lowerseis,'bo')
plot(neqs,upperseis,'bo')
end