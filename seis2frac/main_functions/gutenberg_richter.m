function [n,Bvalue,avalue,btrend,mdl]=gutenberg_richter(mags,bins,Bvalue,avalue,MOC,plot_flg,remove)
%% GUTENBERG_RICHTER
% Script to create GR plots from magnitude data. A- and B-values can be
% provided, or are calculated with a linear fit. If calculating A- and
% B-values from dataset, include a Magnitude of Completion, or it will try
% to fit a linear trend through the low values
dbstop if error
edges=[bins-((bins(2)-bins(1))/2),bins(end)+((bins(2)-bins(1))/2)];

[n]=histcounts(mags,edges);

if isempty('MOC')
    MOC=-Inf;
end

if nargin<7
    remove=0;
end

ordered_orig=sort(mags); % Sort events into ascending order of magnitude
cumfreq_orig=(length(ordered_orig):-1:1); % Assign cumulative frequency to events

% Experimental, from data_density, trying removing the highest magnitude
% bin for extrapolating b-values
ordered=ordered_orig(1:end-remove);
cumfreq=cumfreq_orig(1:end-remove);


ii=find(ordered>=MOC); % Using data only from above Magnitude of Completeness
if size(ii,1)>3 % If there are more than 3 data points, fit linear model
    mdl=fitlm(ordered(ii),log10(cumfreq(ii))); % Linear Model Parameters (decide whether to use them in the next section)
    
    if isempty(avalue) % If avalue was undefined
        if isempty(Bvalue)
            avalue=mdl.Coefficients.Estimate(1); % If bvalue also undefined, use values from linear model
        else % Else check to see if modelled b-value is within 5% of the defined value (assumes input b is correct, and modelled b can be wrong due to too little data)
            if (-mdl.Coefficients.Estimate(2))/Bvalue > 0.95 && (-mdl.Coefficients.Estimate(2))/Bvalue < 1.05
                avalue=mdl.Coefficients.Estimate(1); % If so, use modelled a-value
            else
                avalue=median(log10(cumfreq(ii))'+Bvalue*ordered(ii)); % If not, calculate avalue based of median magnitudes and known bvalue
            end
        end
    end
    
    
    if isempty(Bvalue) % If bvalue hasn't been defined, use value from linear model
        Bvalue=-mdl.Coefficients.Estimate(2);
        fprintf('Calculating b-value from data\n')
    end
    
    
else
    if isempty(Bvalue) % If bvalue hasn't been defined, use value from linear model
        Bvalue=0.85;
        fprintf('Too few data points. B-value hard coded as 0.85\n')
    end
    avalue=log10(mean(cumfreq(ii)))+Bvalue*mean(ordered(ii));
end

btrend=bins*-Bvalue+avalue; % Calculate G-R relation

if plot_flg == 1
    %%
    figure
    figname='G-R Plot';
    set(gcf,'renderer','zbuffer','name',figname); title('Gutenburg-Richter Relation of Relocated Events');
    hold on
    bar(bins,n);
    try
        ylim([0 max(n)/(floor((log10(max(cumfreq))/btrend(1))*10)/10)])
    catch
        fprintf('Error setting Y-axis limits\n')
    end
    yl=ylabel('Frequency');
    xlabel('Magnitude');
    yyaxis right
    plot(ordered,log10(cumfreq),'.'); %ylim([0 btrend(1)])
    if remove > 0
        plot(ordered_orig(end-remove:end),log10(cumfreq_orig(end-remove:end)),'r.');
    end
    plot(bins,btrend,'k--'); ylabel('Log_1_0(Cumulative Frequency)')
    text(max(bins),max(btrend),sprintf('A-value = %.2f',avalue),'HorizontalAlignment','right','VerticalAlignment','bottom')
    text(max(bins),max(btrend),sprintf('B-value = %.2f',Bvalue),'HorizontalAlignment','right','VerticalAlignment','top')
    xline(MOC);
    hold off
    set(yl,'color',[0 0.4470 0.7410])
    
end