function [thickness,fracture_time,ex]=exhumation_time(quadcenters,exhum,seis_depths,seis_ix,QT,plot_fig,fault,bearing,origin,ll)
%% EXHUMATION_TIME
% Using gridded nodal exhumation data, and limits of upper and lower
% seismicity, calculate the amount of time that a rock has to fracture

[ex]=griddata(exhum(:,1),exhum(:,2),exhum(:,3),quadcenters(:,1),quadcenters(:,2));
ex(isnan(ex))=0.1;
[thickness]=griddata(seis_depths(seis_ix,1),seis_depths(seis_ix,2),(seis_depths(seis_ix,3)-seis_depths(seis_ix,4)),quadcenters(:,1),quadcenters(:,2));
fracture_time=(thickness*1e3)./(ex*1e-3);

if plot_fig==1;
    %%
[fault(:,1:2)]=local2llh(fault(:,1:2)',origin)';
for ii=1:QT.BinCount
%     [LL.BinCorners(:,1:3,ii)]=point_rotation(QT.BinCorners(:,1:2,ii),bearing-90,'strike');
    [LL.BinCorners(:,1:2,ii)]=local2llh(QT.BinCorners(:,1:2,ii)',origin)';
    
end
    
    
    
    load C:/Jacko/scripts/ScientificColourMaps6/bamako/bamako.mat; colormap(bamako);
%     bamako=colormap(jet);
    if ~exist('ll')
        [ll(:,1),ll(:,2)]=get_coast(168,173.25,-45,-41.75,[],'f','C:\Jacko\scripts\coastlines\gshhg-bin-2.3.7');
    end
    
  %%  
    figure
    title('Exhumation Rate','FontSize',13)
%     xlabel('Longitude');ylabel('Latitude');
    hold on
    colormap(bamako);colorbar; %caxis([0 25])
    c=colorbar;
    c.Label.String='mm/yr';
    c.Label.FontSize=13;
    caxis([0 8])

    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(LL.BinCorners(:,1,ii),LL.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*(ex(ii)/8)),:))
        end
    end
    plot(ll(:,1),ll(:,2),'k')
                try
        plot(fault([1 2],1),fault([1 2],2),'r');
                end
        xlim([ 168.1365  172.1054]);ylim([ -44.7979  -42.0618])
    %%
    figure
    title('Upper Cut-off')
    xlabel('Longitude');ylabel('Latitude');
    hold on
    colormap(bamako);colorbar; caxis([0 15])
    c=colorbar;

   
    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(LL.BinCorners(:,1,ii),LL.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*seis_depths(ii,3)/-15),:))
        catch
            fill(LL.BinCorners(:,1,ii),LL.BinCorners(:,2,ii),bamako(end,:))
        end
    end
    plot(ll(:,1),ll(:,2),'k')
                try
        plot(fault([1 2],1),fault([1 2],2),'r');
                end
        xlim([ 168.1365  172.1054]);ylim([ -44.7979  -42.0618])
    %%
    figure
    title('Lower Cut-off')
    xlabel('Longitude');ylabel('Latitude');
    hold on
    colormap(bamako);colorbar; caxis([0 25])
    c=colorbar;

    
    try
        plot(fault([1 2],1),fault([1 2],2));
    end
    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(LL.BinCorners(:,1,ii),LL.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*seis_depths(ii,4)/-25),:))
        catch
            if ~isnan(seis_depths(ii,4))
                fill(LL.BinCorners(:,1,ii),LL.BinCorners(:,2,ii),bamako(size(bamako,1),:))
            end
        end
    end

    plot(ll(:,1),ll(:,2),'k')
                try
        plot(fault([1 2],1),fault([1 2],2),'r');
                end
        xlim([ 168.1365  172.1054]);ylim([ -44.7979  -42.0618])
    %%
    figure
    title('Seismogenic Thickness','FontSize',13)
%     xlabel('Longitude');ylabel('Latitude');
    hold on
    colormap(bamako);colorbar; caxis([0 25])
    c=colorbar;
    c.Label.String='km';
c.Label.FontSize=13;

    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(LL.BinCorners(:,1,ii),LL.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*(seis_depths(ii,4)-seis_depths(ii,3))/-25),:))
        catch
            fill(LL.BinCorners(:,1,ii),LL.BinCorners(:,2,ii),bamako(255,:))
        end
    end

    plot(ll(:,1),ll(:,2),'k')
            try
        plot(fault([1 2],1),fault([1 2],2),'r');
            end
        xlim([ 168.1365  172.1054]);ylim([ -44.7979  -42.0618])
    %%
    figure
    title('Time in Seismogenic Zone','FontSize',13)
%     xlabel('Longitude');ylabel('Latitude');
    hold on
    colormap(bamako);colorbar; %caxis([0 25])
    c=colorbar;
    c.Label.String='Million Years';
    caxis([0 5])
    axis equal
    c.Label.FontSize=13;
    for ii=1:QT.BinCount
        try
            fill(LL.BinCorners(:,1,ii),LL.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*(fracture_time(ii)/5e6)),:))
            
        catch
            if ~isnan(fracture_time(ii))
                fill(LL.BinCorners(:,1,ii),LL.BinCorners(:,2,ii),bamako(size(bamako,1),:))
            end
        end
    end

plot(ll(:,1),ll(:,2),'k')
            try
        plot(fault([1 2],1),fault([1 2],2),'r');
    end
    xlim([ 168.1365  172.1054]);ylim([ -44.7979  -42.0618])
end
end