function [C,h]=contour_plot(xy,z,XX,YY,zlims,contour_line_spacing,contour_label_spacing,t,c,fault,coast,bamako)

if contour_label_spacing > range(zlims)
    fprintf('Label Spacing too large. Label all lines\n')
    contour_label_spacing=contour_line_spacing;
end


figure
hold on
oversize2=(z<=zlims(2));
scatter(xy(oversize2,1),xy(oversize2,2),20,z(oversize2),'filled')
scatter(xy(oversize2==0,1),xy(oversize2==0,2),20,'rd','filled')
oversize2=double((z<zlims(2)));
oversize2(oversize2==0)=nan;

if zlims(2)>contour_line_spacing
contour_lines=[zlims(1):contour_line_spacing:floor(zlims(2)*contour_line_spacing)/contour_line_spacing];
else
    contour_lines=[];
end
contour_labels=[zlims(1):contour_label_spacing:max(contour_lines)];

if size(XX,2)~=1 % Contour Gridded Data
    [C,h]=contour(XX(1,:),YY(:,1),reshape(z,size(XX)),contour_lines);
else % Contour Scatttered Data
    tri=delaunay(XX,YY);
    [C,h]=tricontour(tri,XX,YY,z,contour_lines);
end    

try
clabel(C,h,contour_labels);
catch
    fprintf('Contour Label Error\n')
end
colorbar

colormap(bamako)
cb=colorbar;
cb.Label.String=c;
try
caxis([zlims(1) zlims(2)]);
end
plot(fault(1:2,1),fault(1:2,2),'r')

plot(coast(:,1),coast(:,2),'b')
xlim([min(xy(:,1))-25 max(xy(:,1))+25])
ylim([min(xy(:,2))-10 max(xy(:,2))+25])
title(t)