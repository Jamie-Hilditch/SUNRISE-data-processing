function [RSK] = RSKderivedepth(RSK, varargin)

%RSKderivedepth - Calculate depth from pressure.
%
% Syntax:  [RSK] = RSKderivedepth(RSK, [OPTION])
% 
% Calculates depth from pressure and adds the channel metadata in the
% appropriate fields. If the data elements already have a 'depth' channel,
% it is replaced. Uses TEOS-10 toolbox if it is installed
% (http://www.teos-10.org/software.htm#1). Otherwise, it is calculated using
% the Saunders & Fofonoff method.  
%
% Inputs: 
%    [Required] - RSK - Structure containing the logger metadata and data
%
%    [Optional] - latitude - Location of the pressure measurement in
%                       decimal degrees. Default is 45. 
%
% Outputs:
%    RSK - RSK structure containing the depth data
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-07-04

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addParameter(p, 'latitude', 45, @isnumeric);
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
latitude = p.Results.latitude;



RSK = addchannelmetadata(RSK, 'Depth', 'm');
Dcol = getchannelindex(RSK, 'Depth');
[RSKsp, SPcol] = getseapressure(RSK);



castidx = getdataindex(RSK);
for ndx = castidx
    seapressure = RSKsp.data(ndx).values(:, SPcol);
    depth = calculatedepth(seapressure, latitude);
    RSK.data(ndx).values(:,Dcol) = depth;
end



logentry = ['Depth calculated using a latitude of ' num2str(latitude) ' degrees.'];
RSK = RSKappendtolog(RSK, logentry);

end