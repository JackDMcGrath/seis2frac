function [ dmap ] = dataDensity( x, y, width, height, limits, fudge )
%DATADENSITY Get a data density image of data
%   x, y - two vectors of equal length giving scatterplot x, y co-ords
%   width, height - dimensions of the data density plot, in pixels
%   limits - [xmin xmax ymin ymax] - defaults to data max/min
%   fudge - the amount of smear, defaults to size of pixel diagonal
%
% By Malcolm McLean
%
if(nargin == 4)
    limits(1) = min(x);
    limits(2) = max(x);
    limits(3) = min(y);
    limits(4) = max(y);
end
deltax = (limits(2) - limits(1)) / width;
deltay = (limits(4) - limits(3)) / height;
if(nargin < 6)
    fudge = sqrt(deltax^2 + deltay^2);
end
dmap = zeros(height, width);

len = length(x);
% two different vectorisation
if len>height*width
    
    for ii = 0: height - 1
        yi = limits(3) + ii * deltay + deltay/2;
        for jj = 0 : width - 1
            xi = limits(1) + jj * deltax + deltax/2;
            dist2 = (x - xi).^2 + (y - yi).^2;
            dmap(ii+1,jj+1) = sum( 1 ./ ( dist2 + fudge));
        end
    end
    
else
    
    ii = [0:height-1]';
    jj = 0: width-1;
    
    yi = limits(3) + ii * deltay + deltay/2;
    xi = limits(1) + jj * deltax + deltax/2;
    for i = 1:len
        dist2 = (x(i) - xi).^2 + (y(i) - yi).^2;
        ix=(dist2<=10);
%         dmap = dmap + (1 ./ ( dist2 + fudge)) .* ix;
        dmap = dmap + ix;
    end
end

dmap = dmap / (pi*10^2);

end

