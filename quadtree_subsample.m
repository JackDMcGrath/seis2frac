function [QT] = quadtree_subsample(eqs,fault,bearing,QuadminbinCapacity )


Rotation_matrix=rotz(90-bearing);
eqs = (Rotation_matrix\eqs')'; % performing the rotation of the data
fault = (Rotation_matrix\fault')';
[QT]=QuadTree(eqs(:,[1:2]),'minbinCapacity',QuadminbinCapacity,'maxbinCapacity',150,'minSize',10,'grdShape','Square');
%%

%     figure;
%         figname='QuadTree Rotated Subsample';
%
%     set(gcf,'renderer','zbuffer','name',figname); title(figname);
%     hold on
%     boxH = QT.plot;
%     cols = lines(QT.BinCount);
%     doplot = @(p,varargin)plot(p(:,1),p(:,2),varargin{:});
%     for i = 1:QT.BinCount
%         set(boxH(i),'Color',cols(i,:),'LineWidth', 1)
%         doplot(eqs(QT.PointBins==i,:),'.','Color',cols(i,:))
%     end
%     plot(af([1 2],1),af([1 2],2),'LineWidth',3,'Color','k');
%     xlabel('Lon');ylabel('Lat')


% Remove bins that do not contain an event > MOC
[uniqueBins]=unique(QT.PointBins);
QT.BinCount=length(uniqueBins);
QT.BinBoundaries=QT.BinBoundaries(uniqueBins,:);
QT.BinDepths=QT.BinDepths(:,uniqueBins);
QT.BinParents=QT.BinParents(:,uniqueBins);
[in,ix]=ismember(QT.PointBins,uniqueBins); % Indicies of all events in kept bins
QT.Points=QT.Points(in,:);
QT.PointBins=ix(find(ix));


%% Rotate back
Rotation_matrix=rotz(bearing-90);
eqs = (Rotation_matrix\eqs')'; % performing the rotation of the data
fault = (Rotation_matrix\(fault'))';
QT.Points(:,3)=0;
QT.Points = (Rotation_matrix\QT.Points')';
QT.Points = QT.Points(:,[1 2]);
QT.BinCorners = zeros(5,3,QT.BinCount);
for ii=1:QT.BinCount
    binMinMax=QT.BinBoundaries(ii,:);
    QT.BinCorners(:,[1:2],ii)=binMinMax([ 1 2; 3 2 ; 3 4; 1 4; 1 2]);
    QT.BinCorners(:,:,ii)=(Rotation_matrix\QT.BinCorners(:,:,ii)')';
    QT.BinBoundaries(ii,:)=[min(QT.BinCorners(:,[1:2],ii)),max(QT.BinCorners(:,[1:2],ii))];
    
end
QT.BinCorners(:,3,:) = [];


plotfig=1;
if plotfig==1
figure;

if ~exist('ll')
    [lon,lat]=get_coast(168.5,173,-45,-42,[],'f','C:\Jacko\scripts\coastlines\gshhg-bin-2.3.7');
    ll=ll2utm([lon,lat],[169.0844,-43.9019]);
    
end

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
    doplot3(eqs(QT.PointBins==i,:),'.','Color',cols(i,:))
end
plot(fault([1 2],1),fault([1 2],2),'LineWidth',3,'Color','k');
xlabel('Lon');ylabel('Lat');
plot(ll(:,1),ll(:,2))
xlim(xy(:,1));ylim(xy(:,2));
hold off
end
end