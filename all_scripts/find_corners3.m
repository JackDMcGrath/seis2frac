% function [p1,p2,p3,p4]=find_corners2(center,dim,vector)
% Build square around centre point, and define corners
% p1,p2,p3,p4 are TL, TR, BR, BL respectively
%https://www.toppr.com/guides/maths/three-dimensional-geometry/equation-plane-perpendicular-given-vector-passing-through-given-point/
center=pt;
v=[1,0,1];
x=p1(1);
y=p1(2);
z=p1(3);
%% Equation of perpendicular plane passing through centre
q=v(1)*(x-center(1))+v(2)*(y-center(2))+v(3)*(z-center(3))


% end