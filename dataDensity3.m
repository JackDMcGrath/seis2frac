function [ dmap ] = dataDensity3( x, y, z, width, height, depth, radius, limits, fudge )
%DATADENSITY Get a data density image of data
% x, y, z - three vectors of equal length giving scatterplot x, y, z co-ords
% width, height, depth - dimensions of the data density plot, in pixels
% radius - Distance to search for events (so result is in events/km3)
% limits - [xmin xmax ymin ymax zmin zmax] - defaults to data max/min
% fudge - the amount of smear, defaults to size of pixel diagonal
%
% Adapted from Malcolm McLean's dataDensity.m by Jack McGrath
%

if (nargin < 7)
    radius=10;
end

if(nargin < 8) || isempty(limits)
    limits(1) = min(x(:));
    limits(2) = max(x(:));
    limits(3) = min(y(:));
    limits(4) = max(y(:));
    limits(5) = min(z(:));
    limits(6) = max(z(:));
end

deltax = (limits(2) - limits(1)) / width;
deltay = (limits(4) - limits(3)) / height;
deltaz = (limits(6) - limits(5)) / depth;

if(nargin < 9)
    fudge = sqrt(deltax.^2 + deltay.^2 + deltaz.^2);
end

dmap = zeros(height,width,depth);

for ii = 0: height - 1
    yi = limits(3) + ii * deltay + deltay/2;
    
    for jj = 0 : width - 1
        xi = limits(1) + jj * deltax + deltax/2;
        
        for kk = 0 : depth-1
            zi = limits(5) + kk * deltaz + deltaz/2;
            dd = 0;
            dist2 = (x - xi).^2 + (y - yi).^2 + (z - zi).^2;
            dd = (dist2<=radius); %number events within 10k
            dmap(ii+1,jj+1,kk+1) = sum(dd(:))/(4*pi*radius^2); %density (events/km3)
        end
    end
end
end