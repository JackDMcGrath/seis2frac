function [OT,eqs]=octree_subsample(eqs,binCapacity,gridx,plot_figures,MOC,fault,style,gridShape)
%% OCTREE_SUBSAMPLE
% Script to subsample seismicity data
% Requires OcTree.m (edited by JDM)
% https://au.mathworks.com/matlabcentral/fileexchange/40732-octree-partitioning-3d-points-into-spatial-subvolumes
% To use ^ version, just remove 'Division' option from OcTree function
% Requires rotx from Phased Array System Toolbox
% If not installed, rotation matrices can be made manually https://au.mathworks.com/help/phased/ref/rotx.html
%
% Inputs :
%     eqs              : Matrix of seismic data, with Lon, Lat, Depth, Magnitude in columns 1:4
%     binCapacity      : Maximum number of events per bin
%     gridx            : Minimum size of bin
%     plot_figures     : Plot flag
%     MOC (optional)   : Magnitude of Completion (subsampling based on this events above this limit
%     fault (optional) : 4*3 matrix of fault corners (for plotting)
%     style            : OcTree division method ('Normal' or 'Weighted')
%     gridShape        : 'Rectangle' - default subsample or 'Cube' - subsamples as cubes
%
% Outputs
%     OT               : Octree Structure
%     eqs              : Matrix of seismic data inside the bins
%
%     Jack McGrath, University of Leeds, 2021
%     Mar 21 : JDM, Initial commit

rotatestrike=1;
rotatedip=1;
if exist('MOC')~=1
    MOC = -Inf;
end

if exist('fault')~=1 || isempty(fault)
    fault=nan(4,3);
    rotatestrike=0;
    rotatedip=0;
else
    bearing=atand((fault(2,1)-fault(1,1))/(fault(2,2)-fault(1,2)));
    dip=-atand((fault(3,3)-fault(1,3))/sqrt((fault(3,1)-fault(1,1))^2+(fault(3,2)-fault(1,2))^2));
end

if rotatestrike==1 % Rotate data so that grid will be fault parallel
    Rotation_matrix=rotz(90-bearing);
    eqs(:,[1:3]) = (Rotation_matrix\eqs(:,[1:3])')'; % performing the rotation of the data
    fault = (Rotation_matrix\(fault'))';
end

if rotatedip==1 % Rotate data so that grid will be fault parallel
    Rotation_matrix=rotx(dip-90);
    eqs(:,[1:3]) = (Rotation_matrix\eqs(:,[1:3])')'; % performing the rotation of the data
    fault = (Rotation_matrix\(fault'))';
end

OTeqs=eqs(eqs(:,4)>=MOC,:); % Variable to search only for events above MOC

OT = OcTree(OTeqs(:,[1:3]),'binCapacity',binCapacity,'minSize',gridx,'style',style,'grdShape',gridShape); % OcTree sub-sample

% Identify the bins that contain events > MOC
[uniqueBins,~,uBix]=unique(OT.PointBins);

% Add all events to bins
OT.PointBins=OT.query(eqs(:,[1:3])); % Find bins for all data (inc. event < MOC)
OT.Points=eqs(:,[1:3]); % Add locations of all events to variable

% Remove bins that do not contain an event > MOC
OT.BinCount=length(uniqueBins);
OT.BinBoundaries=OT.BinBoundaries(uniqueBins,:);
OT.BinDepths=OT.BinDepths(:,uniqueBins);
OT.BinParents=OT.BinParents(:,uniqueBins);
[in,ix]=ismember(OT.PointBins,uniqueBins); % Indicies of all events in kept bins
OT.Points=OT.Points(in,:);
OT.PointBins=ix(find(ix));
eqs=eqs(in,:);

%%
% if plot_figures==1
%     figure;
%     figname='OcTree Subsample rotated';
%     set(gcf,'renderer','zbuffer','name',figname); title(figname);
%     hold on
%     boxH = OT.plot;
%     cols = lines(OT.BinCount);
%     doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
%     for i = 1:OT.BinCount
%         set(boxH(i),'Color',cols(i,:),'LineWidth', 1)
%         doplot3(eqs(OT.PointBins==i,:),'.','Color',cols(i,:))
%     end
%     plot3(fault([1 2 4 3 1],1),fault([1 2 4 3 1],2),fault([1 2 4 3 1],3),'LineWidth',3,'Color','k');
%     xlabel('Parallel');ylabel('Perpendicular');zlabel('Depth');pbaspect([1 1 1]);
%     axis image, view([-90,0])
% end

if rotatedip==1 % Rotate data back into unrotated form
    Rotation_matrix=rotx(90-dip);
    eqs(:,[1:3]) = (Rotation_matrix\eqs(:,[1:3])')'; % performing the rotation of the data
    fault = (Rotation_matrix\(fault'))';
    OT.Points = (Rotation_matrix\OT.Points')';
    for ii=1:OT.BinCount
        binMinMax=OT.BinBoundaries(ii,:);
        OT.BinCorners(:,:,ii)=binMinMax([ 1 2 3; 4 2 3; 4 5 3; 1 5 3;
            1 2 6; 4 2 6; 4 5 6; 1 5 6]);
        OT.BinCorners(:,:,ii)=(Rotation_matrix\OT.BinCorners(:,:,ii)')';
        OT.BinBoundaries(ii,:)=[min(OT.BinCorners(:,:,ii)),max(OT.BinCorners(:,:,ii))];
    end
    
    if rotatestrike==1 % Rotate data back into unrotated form
        Rotation_matrix=rotz(bearing-90);
        eqs(:,[1:3]) = (Rotation_matrix\eqs(:,[1:3])')'; % performing the rotation of the data
        fault = (Rotation_matrix\(fault'))';
        OT.Points = (Rotation_matrix\OT.Points')';
        
        for ii=1:OT.BinCount
            OT.BinCorners(:,:,ii)=(Rotation_matrix\OT.BinCorners(:,:,ii)')';
            OT.BinBoundaries(ii,:)=[min(OT.BinCorners(:,:,ii)),max(OT.BinCorners(:,:,ii))];
        end
        
    end
    
elseif rotatestrike==1 % Rotate data back into unrotated form
    Rotation_matrix=rotz(bearing-90);
    eqs(:,[1:3]) = (Rotation_matrix\eqs(:,[1:3])')'; % performing the rotation of the data
    fault = (Rotation_matrix\(fault'))';
    OT.Points = (Rotation_matrix\OT.Points')';
    for ii=1:OT.BinCount
        binMinMax=OT.BinBoundaries(ii,:);
        OT.BinCorners(:,:,ii)=binMinMax([ 1 2 3; 4 2 3; 4 5 3; 1 5 3;
            1 2 6; 4 2 6; 4 5 6; 1 5 6]);
        OT.BinCorners(:,:,ii)=(Rotation_matrix\OT.BinCorners(:,:,ii)')';
        OT.BinBoundaries(ii,:)=[min(OT.BinCorners(:,:,ii)),max(OT.BinCorners(:,:,ii))];
    end
    
end

if plot_figures == 1
    
    figure;
    if rotatestrike==1 || rotatedip==1
        figname='OcTree Subsample unrotated';
    else
        figname='OcTree Subsample';
    end
    set(gcf,'renderer','zbuffer','name',figname); title(figname);
    hold on
    boxH = OT.plot;
    cols = lines(OT.BinCount);
    doplot3 = @(p,varargin)plot3(p(:,1),p(:,2),p(:,3),varargin{:});
    for i = 1:OT.BinCount
        set(boxH(i),'Color',cols(i,:),'LineWidth', 1)
        doplot3(eqs(OT.PointBins==i,:),'.','Color',cols(i,:))
    end
    plot3(fault([1 2 4 3 1],1),fault([1 2 4 3 1],2),fault([1 2 4 3 1],3),'LineWidth',3,'Color','k');
    xlabel('Lon');ylabel('Lat');zlabel('Depth');pbaspect([1 1 1]);
    axis image,  view([-36,25])
    
end

end
