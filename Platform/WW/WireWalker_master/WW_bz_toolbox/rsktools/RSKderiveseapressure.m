function [RSK] = RSKderiveseapressure(RSK, varargin)

%RSKderiveseapressure - Calculate sea pressure.
%
% Syntax:  [RSK] = RSKderiveseapressure(RSK, [OPTIONS])
% 
% Derives sea pressure and fills all of data's elements and channel
% metadata. If sea pressure already exists, it recalculates it and
% overwrites that data column.  
%
% Inputs: 
%    [Required] - RSK - Structure containing the logger metadata and data
%
%    [Optional] - patm - Atmospheric Pressure. Default is value stored in
%                       parameters table or 10.1325 dbar if unavailable. 
%
% Outputs:
%    RSK - Structure containing the sea pressure data.
%
% See also: getseapressure, RSKplotprofiles.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-07-04

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addParameter(p, 'patm', [], @isnumeric);
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
patm = p.Results.patm;


try
    Pcol = getchannelindex(RSK, 'Pressure');
catch
    Pcol = getchannelindex(RSK, 'BPR pressure');
end

if isempty(patm)
    patm = getatmosphericpressure(RSK);
end



RSK = addchannelmetadata(RSK, 'Sea Pressure', 'dbar');
SPcol = getchannelindex(RSK, 'Sea Pressure');



castidx = getdataindex(RSK);
for ndx = castidx
    seapressure = RSK.data(ndx).values(:, Pcol)- patm;
    RSK.data(ndx).values(:,SPcol) = seapressure;
end



logentry = ['Sea pressure calculated using an atmospheric pressure of ' num2str(patm) ' dbar.'];
RSK = RSKappendtolog(RSK, logentry);

end


