function [gps]=FPGPS(gps,bearing,px)
% Convert GPS velocities to Fault Parallel (bearing measured anticlockwise
% from east)

gps.vx=[gps.ve*cosd(bearing)+gps.vn*sind(bearing)];
gps.vy=[gps.ve*sind(bearing)-gps.vn*cosd(bearing)];

% Have to subtract the perpendicular (imagine ve +1 and vn +1 on a 45 degree
% fault- the fault perps cancel each other out). Validated by checking that
% the combined EN velocity matches combined XY