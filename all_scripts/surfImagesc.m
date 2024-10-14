function g =  surfImagesc(x,y,z,varargin)
% This function plots the 2d imagesc plot along with 3-D surf plot
% x = x coordinate (vector)
% y = y coordinate (vector)
% z = Z coordinate (The value to be plotted)
% cmType = type of predefined color map (jet,parula, hot, gray etc.)--
% refer https://in.mathworks.com/help/matlab/ref/colormap.html 
% clim=[lowerBound higherBound]-- Limit for the colormap
% vAngle=[Azimuth Elevation]-- viewing angle
% offset = movement of image along z-axis
if nargin < 3
    error('Not enoungh input arguments. This function expects atleast 3 input arguments')
end
   defaultcmType = jet(50);
   defaultCLim = [min(min(z)) max(max(z))];
   defaultvAngle = [-45 30];
   defaultOffset = min(min(z));
   p = inputParser;
   
   validVec = @(x) isnumeric(x) && isvector(x);
   validMAT = @(x) isnumeric(x) && ismatrix(x);
   validScalar = @(x) isscalar(x);
   addRequired(p,'x',validVec);
   addRequired(p,'y',validVec);
   addRequired(p,'z',validMAT);
   addParameter(p,'cmType',defaultcmType ,validMAT);
   addParameter(p,'CLim',defaultCLim ,validMAT);
   addParameter(p,'vAngle',defaultvAngle ,validVec);
   addParameter(p,'Offset',defaultOffset ,validScalar);
    
   parse(p,x,y,z,varargin{:});
   surf(p.Results.x,p.Results.y,p.Results.z);hold on
   
   g = hgtransform('Matrix',makehgtform('translate',[0 0 p.Results.Offset]));
   imagesc(g,p.Results.x,p.Results.y,p.Results.z,p.Results.CLim);
   colormap(p.Results.cmType);
    xlabel('x');
    ylabel('y');
    colorbar;
    view(p.Results.vAngle);
    hold off;
end