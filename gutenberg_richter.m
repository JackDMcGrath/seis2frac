function [n,bvalue,avalue,btrend,mdl]=gutenberg_richter(mags,bins,bvalue,avalue,MOC,plot_flg)
%% GUTENBERG_RICHTER
% Script to create GR plots from magnitude data. A- and B-values can be
% provided, or are calculated with a linear fit. If calculating A- and
% B-values from dataset, include a Magnitude of Completion, or it will try
% to fit a linear trend through the low values

edges=[bins-((bins(2)-bins(1))/2),bins(end)+((bins(2)-bins(1))/2)];
if plot_flg==1;
    plot_flg;end
[n]=histcounts(mags,edges);

if isempty('MOC')
    MOC=-Inf;
end

ordered=sort(mags); % Sort events into ascending order of magnitude
cumfreq=(length(ordered):-1:1); % Assign cumulative frequency to events

ii=find(ordered>=MOC); % Using data only from above Magnitude of Completeness
if size(ii,1)>3 % If there are more than 3 data points, fit linear model
    mdl=fitlm(ordered(ii),log10(cumfreq(ii))); % Linear Model Parameters
    
    if isempty(avalue) % If avalue was undefined
        if isempty(bvalue)
            avalue=mdl.Coefficients.Estimate(1); % If bvalue also undefined, use values from linear model
        else % Else check to see if modelled b-value is within 5% of the defined value (assumes input b is correct, and modelled b can be wrong due to too little data)
            if (-mdl.Coefficients.Estimate(2))/bvalue > 0.95 && (-mdl.Coefficients.Estimate(2))/bvalue < 1.05
                avalue=mdl.Coefficients.Estimate(1); % If so, use modelled a-value
            else
                avalue=log10(mean(cumfreq(ii)))+bvalue*mean(ordered(ii)); % If not, calculate avalue based of mean magnitudes and known bvalue
            end
        end
    end
    
    
    if isempty(bvalue) % If bvalue hasn't been defined, use value from linear model
        bvalue=-mdl.Coefficients.Estimate(2);
    end
    
    
else
    avalue=log10(mean(cumfreq(ii)))+bvalue*mean(ordered(ii));
end

btrend=bins*-bvalue+avalue; % Calculate G-R relation

if plot_flg == 1
    %%
    figure
    figname='G-R Plot';
    set(gcf,'renderer','zbuffer','name',figname); title('Gutenburg-Richter Relation of Relocated Events');
    hold on
    bar(bins,n); ylim([0 max(n)/(floor((log10(max(cumfreq))/btrend(1))*10)/10)])
    xlabel('Magnitude'); ylabel('Frequency')
    yyaxis right
    plot(ordered,log10(cumfreq),'.'); %ylim([0 btrend(1)])
    plot(bins,btrend,'k--'); ylabel('Log_1_0(Cumulative Frequency)')
    xline(MOC);
    hold off
end