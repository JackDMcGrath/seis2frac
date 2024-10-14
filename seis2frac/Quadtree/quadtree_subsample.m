function [QT,dropped,ll] = quadtree_subsample(EQS,FAULT,bearing,dip,QuadminbinCapacity,r_flg,ll,shrink,plotfig,grdshp,maxsize)
%% QUADTREE_SUBSAMPLE
% Function to carry out quadtree subsampling of point data, with the option
% of rotating the data to a given strike and dip in order to optimise the
% gridsearch.
%
% Jack McGrath, Apr 2021

% plotfig=1;


if ~exist('ll') && plotfig==1
    [lon,lat]=get_coast(168.5,173,-45,-42,[],'f','C:\Jacko\scripts\coastlines\gshhg-bin-2.3.7');
    ll=ll2utm([lon,lat],[169.0844,-43.9019]);
    ll(:,3)=0;
end

if r_flg.norm2strike==1
    [eqs]=point_rotation(EQS(:,1:3),90-bearing,'strike');
    [fault]=point_rotation(FAULT,90-bearing,'strike');
else
    eqs=EQS(:,1:3);
    fault=FAULT;
end

if plotfig~=0
    figure;
    figname='Seismicity Distribution';
    set(gcf,'renderer','zbuffer','name',figname); title(figname);
    hold on
    scatter3(eqs(:,1),eqs(:,2),eqs(:,3),5,EQS(:,4),'filled');
    colorbar
    load C:/Jacko/scripts/ScientificColourMaps6/bamako/bamako.mat;colormap(bamako);
    try
        plot3(fault([1 2 4 3 1],1),fault([1 2 4 3 1],2),fault([1 2 4 3 1],3),'r');
    end
    if r_flg.norm2strike==1
        xlabel('Along Strike Distance (km)');ylabel('Fault Perpendicular Distance');zlabel('Depth');
        mtcook=point_rotation([117.4818,56.1926],90-bearing,'strike');
        plot3(mtcook(1),mtcook(2),3.5,'kp','MarkerFaceColor','g','MarkerEdgeColor','w','MarkerSize',12)
        [ll_rot]=point_rotation(ll,90-bearing,'strike');
        plot(ll_rot(:,1),ll_rot(:,2))
    else
        xlabel('Lon');ylabel('Lat');zlabel('Depth');
        plot3(117.4818,56.1926,3.5,'kp','MarkerFaceColor','g','MarkerEdgeColor','w','MarkerSize',12)
        plot(ll(:,1),ll(:,2))
        ll_rot=ll;
    end
    axis equal;view([0 0])
    hold off
end


if r_flg.norm2dip==1 % Rotate data so that grid will be fault parallel
    if r_flg.norm2strike==0
        fprintf('WARNING: Rotating Dip without first rotating to strike. Results could be funky! \n')
    end
    [eqs]=point_rotation(eqs,dip-90,'dip');
    [fault]=point_rotation(fault,dip-90,'dip');
    
end

if r_flg.norm2strike==0 && r_flg.norm2dip==0 % If no rotation is to be carried out
    eqs=EQS;
    fault=FAULT;
end

[QT]=QuadTree(eqs(:,[1:2]),'minbinCapacity',QuadminbinCapacity,'maxbinCapacity',3*QuadminbinCapacity,'minSize',10,'grdShape',grdshp,'maxSize',maxsize);

if strcmpi(shrink,'shrink')
    QT.shrink;
end
%%
if plotfig==1
    figure;
    figname='QuadTree Rotated Subsample';
    
    set(gcf,'renderer','zbuffer','name',figname); title(figname);
    hold on
    boxH = QT.plot;
    cols = lines(QT.BinCount);
    doplot = @(p,varargin)plot(p(:,1),p(:,2),varargin{:});
    doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
    for i = 1:QT.BinCount
        set(boxH(i),'Color',cols(i,:),'LineWidth', 1)
        try
            doplot3(eqs(QT.PointBins==i,:),'.','Color',cols(i,:))
        catch
            doplot(eqs(QT.PointBins==i,:),'.','Color',cols(i,:))
        end
    end
    plot3(fault([1 2 4 3 1],1),fault([1 2 4 3 1],2),fault([1 2 4 3 1],3),'LineWidth',3,'Color','k');
    xlabel('Lon');ylabel('Lat')
end

% Remove bins that do not contain enough events
uniqueBins=[];
for ii=1:QT.BinCount
    count=sum(QT.PointBins==ii);
    if count>=QuadminbinCapacity
        uniqueBins=[uniqueBins;ii];
    end
end
% [uniqueBins]=unique(QT.PointBins); % Find all bins with an event
QT.BinCount=length(uniqueBins); % Number of Bins to keep
QT.BinBoundaries=QT.BinBoundaries(uniqueBins,:); % Boundaries of said bins
QT.BinDepths=QT.BinDepths(:,uniqueBins);
QT.BinParents=QT.BinParents(:,uniqueBins);
QT.Points(:,3)=eqs(:,3);
[in,ix]=ismember(QT.PointBins,uniqueBins); % Indicies of all events in kept bins
dropped=find(in==0);
QT.Points=QT.Points(in,:);
QT.PointBins=ix(find(ix));
QT.BinCorners = zeros(5,3,QT.BinCount);

%% Rotate back
if r_flg.dip2norm==1
    [QT.Points]=point_rotation(QT.Points,90-dip,'dip');
end

if r_flg.strike2norm==1
    [QT.Points]=point_rotation(QT.Points,bearing-90,'strike');
    %     [QT.Points]=EQS(in,:);
    for ii=1:QT.BinCount
        binMinMax=QT.BinBoundaries(ii,:);
        QT.BinCorners(:,[1:2],ii)=binMinMax([ 1 2; 3 2 ; 3 4; 1 4; 1 2]);
        [QT.BinCorners(:,:,ii)]=point_rotation(QT.BinCorners(:,:,ii),bearing-90,'strike');
        QT.BinBoundaries(ii,:)=[min(QT.BinCorners(:,[1:2],ii)),max(QT.BinCorners(:,[1:2],ii))];
    end
else
    for ii=1:QT.BinCount
        binMinMax=QT.BinBoundaries(ii,:);
        QT.BinCorners(:,[1:2],ii)=binMinMax([ 1 2; 3 2 ; 3 4; 1 4; 1 2]);
    end
end

if plotfig==1
    figure;
    
    
    xy=[168.75,-44.5; 171.75,-42.5];
    xy=ll2utm(xy,[169.0844,-43.9019]);
    
    figname='QuadTree Subsample';
    
    set(gcf,'renderer','zbuffer','name',figname); title(figname);
    hold on
    boxH = QT.plot;
    cols = lines(QT.BinCount);
    doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
    for i = 1:QT.BinCount
        set(boxH(i),'Color',cols(i,:),'LineWidth', 1)
        doplot3(QT.Points(QT.PointBins==i,:),'.','Color',cols(i,:))
    end
    if r_flg.strike2norm==0
        xlabel('Along Strike Distance (km)');ylabel('Fault Perpendicular Distance');zlabel('Depth');
        mtcook=point_rotation([117.4818,56.1926],90-bearing,'strike');
        plot3(mtcook(1),mtcook(2),3.5,'kp','MarkerFaceColor','g','MarkerEdgeColor','w','MarkerSize',12)
        plot(ll_rot(:,1),ll_rot(:,2))
        plot3(fault([1 2 4 3 1],1),fault([1 2 4 3 1],2),fault([1 2 4 3 1],3),'LineWidth',3,'Color','k');
    else
        xlabel('Lon');ylabel('Lat');zlabel('Depth');
        plot3(117.4818,56.1926,3.5,'kp','MarkerFaceColor','g','MarkerEdgeColor','w','MarkerSize',12)
        plot(ll(:,1),ll(:,2))
        plot3(FAULT([1 2 4 3 1],1),FAULT([1 2 4 3 1],2),FAULT([1 2 4 3 1],3),'LineWidth',3,'Color','k');
    end
    xlim(xy(:,1));ylim(xy(:,2));
    hold off
    axis equal
    xlim([floor(min(QT.Points(:,1))/50)*50 ceil(max(QT.Points(:,1))/50)*50]);ylim([floor(min(QT.Points(:,2))/50)*50 ceil(max(QT.Points(:,2))/50)*50])
end
end