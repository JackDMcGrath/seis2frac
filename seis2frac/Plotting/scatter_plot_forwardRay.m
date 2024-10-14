function contour_plot_forwardRay(xy,data,zlims,t,c,fault,coast,bamako)

figure
hold on
oversize=[];

% Plot Boxes
for ii=1:size(data,1)
    if isnan(data(ii))
        scatter(xy(ii,1),xy(ii,2),1500,0,'c','filled')
%         text(xy(ii,1),xy(ii,2),sprintf('%.0f',data(ii)),'Color','g','HorizontalAlignment','Center')
    elseif data(ii)>zlims(2)
        scatter(xy(ii,1),xy(ii,2),1500,zlims(2),'d','filled')
%         text(xy(ii,1),xy(ii,2),sprintf('%.0f',data(ii)),'Color','g','HorizontalAlignment','Center')
    else
        scatter(xy(ii,1),xy(ii,2),1500,data(ii),'s','filled')
%         text(xy(ii,1),xy(ii,2),sprintf('%.0f',data(ii)),'Color','g','HorizontalAlignment','Center')
    end
end

% Plot Labels
for ii=1:size(data,1)
    if isnan(data(ii))
%         scatter(xy(ii,1),xy(ii,2),1500,0,'c','filled')
        text(xy(ii,1),xy(ii,2),sprintf('%.0f',data(ii)),'Color','g','HorizontalAlignment','Center')
    elseif data(ii)>zlims(2)
%         scatter(xy(ii,1),xy(ii,2),1500,zlims(2),'d','filled')
        text(xy(ii,1),xy(ii,2),sprintf('%.0f',data(ii)),'Color','g','HorizontalAlignment','Center')
    else
%         scatter(xy(ii,1),xy(ii,2),1500,data(ii),'s','filled')
        text(xy(ii,1),xy(ii,2),sprintf('%.0f',data(ii)),'Color','g','HorizontalAlignment','Center')
    end
end

colormap(bamako)
cb=colorbar;
cb.Label.String=c;
try
caxis([zlims(1) zlims(2)]);
end
plot(fault(1:2,1),fault(1:2,2),'r')

plot(coast(:,1),coast(:,2),'k')
xlim([min(xy(:,1))-25 max(xy(:,1))+25])
ylim([min(xy(:,2))-10 max(xy(:,2))+25])
title(t)