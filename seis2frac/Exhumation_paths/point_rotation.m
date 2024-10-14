function [rpts]=point_rotation(pts,angle,axis)
%% POINT_ROTATION
% Function to rotate a set of points around a given:
%     strike (Rotate around z axis)
%     dip (Rotate around x axis [must have rotated to strike first)
%     rake (Rotate around y axis)

% Jack McGrath, Apr 2021

if size(pts,2)~=3 % Quick check to make sure 2D data is given a 3rd dimension
    pts(:,3)=0;
end


if strcmpi(axis,'strike')
    Rotation_matrix=rotz(angle);
    rpts = (Rotation_matrix\pts')'; % performing the rotation of the data
end

if strcmpi(axis,'dip')
    Rotation_matrix=rotx(angle);
    rpts = (Rotation_matrix\pts')'; % performing the rotation of the data
end

if strcmpi(axis,'rake')
    Rotation_matrix=roty(angle);
    rpts = (Rotation_matrix\pts')'; % performing the rotation of the data
end

end