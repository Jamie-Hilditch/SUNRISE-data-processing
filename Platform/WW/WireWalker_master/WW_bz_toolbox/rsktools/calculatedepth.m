function depth = calculatedepth(pressure, latitude)

%CALCULATEDEPTH - Calculate depth from pressure.
%
% Syntax:  [depth] = CALCULATEDEPTH(pressure, latitude)
% 
% Calculates depth using pressure data and latitude. Uses TEOS-10 toolbox
% if it is installed. The toolbox can be found at 
% http://www.teos-10.org/software.htm#1. Otherwise, it is calculated using
% the Saunders & Fofonoff method.  
% 
% Inputs:
%    pressure - Vector of pressure values in dbar
%
%    latitude - Location of the pressure measurement in decimal degrees
%               north. 
%
% Outputs:
%    depth - Vector containing depth in meters.
%
% Example: 
%    depth = CALCULATEDEPTH(pressure, 52)
%
% See also: RSKderivedepth.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-20

hasTEOS = ~isempty(which('gsw_z_from_p'));

if hasTEOS
    depth = -gsw_z_from_p(pressure, latitude);  
    
else
    x = (sin(latitude/57.29578)).^2;
    gr = 9.780318*(1.0 + (5.2788e-3 + 2.36e-5*x).*x) + 1.092e-6.*pressure;
    depth = (((-1.82e-15*pressure + 2.279e-10).*pressure - 2.2512e-5).*pressure + 9.72659).*pressure;
    depth = depth./gr;
end

end