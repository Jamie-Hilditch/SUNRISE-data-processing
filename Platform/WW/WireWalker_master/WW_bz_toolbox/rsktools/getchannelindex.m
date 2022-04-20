function channelIdx = getchannelindex(RSK, channel)

%GETCHANNELINDEX - Return index of channels.
%
% Syntax:  [channelIdx] = GETCHANNELINDEX(RSK, channel)
% 
% Finds the channel index in the RSK of the channel longNames given. If the
% channel is not in the RSK, it returns an error.
%
% Inputs:
%   RSK - RSK structure
%
%   channel - LongName as written in RSK.channels.
%
% Outputs:
%    channelIdx - Array containing the index of channels.
%
% See also: RSKplotdata, RSKsmooth, RSKderivedepth.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-07-04

if any(strcmpi(channel, {RSK.channels.longName}));
    chanCol = find(strcmpi(channel, {RSK.channels.longName}));
    channelIdx = chanCol(1);
else
    error(['There is no ' channel ' channel in this file.']);
end

end