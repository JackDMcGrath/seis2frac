function [m]=rot3d(vector)
%Function that will generate a rotation matrix to rotate a plane to
%perpendicular to a given vector


zangle=atand(vector(2)/vector(1)); % Rotation angle in XY plane
yangle=atand(vector(1)/vector(3)); % Rotation angle in XZ plane
xangle=atand(vector(2)/vector(3)); % Rotation angle in YZ plane

%% Checks in the event of no rotation needed in a plane
if isnan(zangle)
    zangle=0;
end

if isnan(yangle)
    yangle=0;
end

if isnan(xangle)
    xangle=0;
end

m=rotz(-zangle)*roty(-yangle)*rotx(xangle); % Create 3-D rotation matrix

end