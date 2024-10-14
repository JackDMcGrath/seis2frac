function [thickness,fracture_time,ex]=exhumation_time_surface(octcenters,quadcenters,exhum,seis_depths,seis_ix,OT,QT,plot_fig,fault,origin,bearing,dip)

% octcenters= bin_center;
% quadcenters=quad_center;
% exhum=exhum;
% plot_fig=1;
% fault=af;


%% QT = output from Octree or Quadtree

[ex]=griddata(exhum(:,1),exhum(:,2),exhum(:,3),octcenters(:,1),octcenters(:,2));
[qex]=griddata(exhum(:,1),exhum(:,2),exhum(:,3),quadcenters(:,1),quadcenters(:,2));
ex(isnan(ex))=0.1;
qex(isnan(ex))=0.1;

thickness=nan(OT.BinCount,1);
for ii=1:OT.BinCount
[thickness(ii)]=depth_surf(octcenters(ii,:),seis_depths,bearing,dip,ii);
end

% [thickness]=griddata(seis_depths(seis_ix,1),seis_depths(seis_ix,2),(seis_depths(seis_ix,3)-seis_depths(seis_ix,4)),octcenters(:,1),octcenters(:,2));
fracture_time=(thickness*1e3)./(ex*1e-3);
[qthickness]=griddata(seis_depths(seis_ix,1),seis_depths(seis_ix,2),(seis_depths(seis_ix,3)-seis_depths(seis_ix,4)),quadcenters(:,1),quadcenters(:,2));
qfracture_time=(qthickness*1e3)./(qex*1e-3);

if plot_fig==1
    %%
    load C:/Jacko/scripts/ScientificColourMaps6/bamako/bamako.mat
    % bamako=colormap(jet);
    if ~exist('ll')
        [lon,lat]=get_coast(169,173,-45,-42,[],'f','C:\Jacko\scripts\coastlines\gshhg-bin-2.3.7');
        ll=ll2utm([lon,lat],origin);
    end
    
    %%
    figure
    title('Oct Exhumation Rate')
    xlabel('Lon');ylabel('Lat');
    hold on
    plot(ll(:,1),ll(:,2))
    colormap(bamako);colorbar; %caxis([0 25])
    c=colorbar;
    caxis([0 10])
    
    try
        plot(fault([1 2],1),fault([1 2],2));
    end
    axis equal
    for ii=1:OT.BinCount
        try
            fill(OT.BinCorners([8 7 2 1],1,ii),OT.BinCorners([8 7 2 1],2,ii),bamako(ceil(size(bamako,1)*(ex(ii)/10)),:))
        end
    end
    
        figure
    title('Quad Exhumation Rate')
    xlabel('Lon');ylabel('Lat');
    hold on
    plot(ll(:,1),ll(:,2))
    colormap(bamako);colorbar; %caxis([0 25])
    c=colorbar;
    caxis([0 10])
    
    try
        plot(fault([1 2],1),fault([1 2],2));
    end
    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*(qex(ii)/10)),:))
        end
    end
    %%
    figure
    title('Upper Cut-off')
    xlabel('Lon');ylabel('Lat');
    hold on
    plot(ll(:,1),ll(:,2))
    colormap(bamako);colorbar; caxis([0 25])
    c=colorbar;
    
    try
        plot(fault([1 2],1),fault([1 2],2));
    end
    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*seis_depths(ii,3)/-25),:))
        end
    end
    %%
    figure
    title('Lower Cut-off')
    xlabel('Lon');ylabel('Lat');
    hold on
    plot(ll(:,1),ll(:,2))
    colormap(bamako);colorbar; caxis([0 25])
    c=colorbar;
    
    try
        plot(fault([1 2],1),fault([1 2],2));
    end
    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*seis_depths(ii,4)/-25),:))
        catch
            if ~isnan(seis_depths(ii,4))
                fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(size(bamako,1),:))
            end
        end
    end
    %%
    figure
    title('Thickness')
    xlabel('Lon');ylabel('Lat');
    hold on
    plot(ll(:,1),ll(:,2))
    colormap(bamako);colorbar; caxis([0 25])
    c=colorbar;
    
    try
        plot(fault([1 2],1),fault([1 2],2));
    end
    axis equal
    
    for ii=1:QT.BinCount
        try
            fill(QT.BinCorners(:,1,ii),QT.BinCorners(:,2,ii),bamako(ceil(size(bamako,1)*(seis_depths(ii,4)-seis_depths(ii,3))/-25),:))
        end
    end
    
        %%
    figure
    title('Oct Thickness')
    xlabel('Lon');ylabel('Lat');
    hold on
    plot(ll(:,1),ll(:,2))
    colormap(bamako);colorbar; %caxis([0 25])
    c=colorbar;
    caxis([0 25])
    
    try
        plot(fault([1 2],1),fault([1 2],2));
    end
    axis equal
    
    for ii=1:OT.BinCount
        try
            fill(OT.BinCorners([8 7 2 1],1,ii),OT.BinCorners([8 7 2 1],2,ii),bamako(ceil(size(bamako,1)*(thickness(ii)/25)),:))
        catch
            if ~isnan(fracture_time(ii))
                fill(OT.BinCorners([8 7 2 1],1,ii),OT.BinCorners([8 7 2 1],2,ii),bamako(size(bamako,1),:))
            end
        end
    end
    
    %%
    figure
    title('Fracture Time')
    xlabel('Lon');ylabel('Lat');
    hold on
    plot(ll(:,1),ll(:,2))
    colormap(bamako);colorbar; %caxis([0 25])
    c=colorbar;
    caxis([0 5e6])
    
    try
        plot(fault([1 2],1),fault([1 2],2));
    end
    axis equal
    
    for ii=1:OT.BinCount
        try
            fill(OT.BinCorners([8 7 2 1],1,ii),OT.BinCorners([8 7 2 1],2,ii),bamako(ceil(size(bamako,1)*(fracture_time(ii)/5e6)),:))
        catch
            if ~isnan(fracture_time(ii))
                fill(OT.BinCorners([8 7 2 1],1,ii),OT.BinCorners([8 7 2 1],2,ii),bamako(size(bamako,1),:))
            end
        end
    end

end
end