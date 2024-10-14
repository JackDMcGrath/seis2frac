function [p1,p2,p3,p4]=find_corners(center,width)
% Build square around centre point, and define corners
% p1,p2,p3,p4 are TL, TR, BR, BL respectively

dim=width/2;

p1=center+[-dim, dim 0];
p2=center+[ dim, dim 0];
p3=center+[ dim,-dim 0];
p4=center+[-dim,-dim 0];

end