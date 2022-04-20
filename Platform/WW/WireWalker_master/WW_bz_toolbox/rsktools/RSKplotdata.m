function handles = RSKplotdata(RSK, varargin)

%RSKplotdata - Plot summaries of logger data.
%
% Syntax:  [handles] = RSKplotdata(RSK, [OPTIONS])
% 
% Generates a plot, similar to the thumbnail plot, only using the full
% 'data' that you read in. If data field has many elements, many casts, the
% default is to use the first data element; RSK.data(1).
% 
% Inputs:
%    [Required] - RSK - Structure containing the logger metadata and data.
%
%    [Optional] - channel - Longname of channel to plot, can be multiple in
%                       a cell, if no value is given it will plot all
%                       channels.
%
%                 profile - Profile number. Default is 1.
% 
%                 direction - 'up' for upcast, 'down' for downcast. Default
%                       is the first string in RSK.profiles.order; the
%                       first cast.
%
% Output:
%     handles - Line object of the plot.
%
% Example: 
%    RSK = RSKopen('sample.rsk');   
%    RSK = RSKreaddata(RSK);  
%    RSKplotdata(RSK);
%    -OR-
%    handles = RSKplotdata(RSK, 'channel', {'Temperature', 'Conductivity'})
%
% See also: RSKplotprofiles, RSKplotburstdata.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-19

validDirections = {'down', 'up'};
checkDirection = @(x) any(validatestring(x,validDirections));

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addParameter(p, 'channel', 'all');
addParameter(p, 'profile', 1);
addParameter(p, 'direction', [], checkDirection);
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
channel = p.Results.channel;
profile = p.Results.profile;
direction = p.Results.direction;



if ~isfield(RSK,'data')
    disp('You must read a section of data in first!');
    disp('Use RSKreaddata...')
    handles = [];
    return
end


castidx = getdataindex(RSK, profile, direction);
if isfield(RSK, 'profiles') && isfield(RSK.profiles, 'order') && any(strcmp(p.UsingDefaults, 'direction'))
    direction = RSK.profiles.order{1};
    castidx = getdataindex(RSK, profile, direction);
end
if size(castidx,2) ~= 1 
    disp('RSKplotdata can only plot one cast, use RSKplotprofiles...');
    return
end
    


chanCol = [];
channels = cellchannelnames(RSK, channel);
for chan = channels
    chanCol = [chanCol getchannelindex(RSK, chan{1})];
end

handles = channelsubplots(RSK, 'data', 'chanCol', chanCol, 'castidx', castidx);

end


