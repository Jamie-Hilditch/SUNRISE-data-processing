function handles = RSKplotdownsample(RSK, varargin)

% RSKplotdownsample - Plot summaries of logger data downsample.
%
% Syntax:  [handles] = RSKplotdownsample(RSK, [OPTIONS])
% 
% Generates a summary plot of the downsample data in the RSK structure.
% 
% Inputs:
%    [Required] - RSK - Structure containing the logger metadata and
%                       downsample.
%
%    [Optional] - channel - Longname of channel to plots, can be multiple
%                           in a cell, if no value is given it will plot
%                           all channels. 
%
% Output:
%    handles - Line object of the plot.
%
% Example: 
%    RSK = RSKopen('sample.rsk');  
%    RSKplotdownsample(RSK);  
%
% See also: RSKopen, RSKplotdata, RSKplotburstdata.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2018-01-24

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addParameter(p, 'channel', 'all');
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
channel = p.Results.channel;



field = 'downsample';
if ~isfield(RSK,field)
    disp('You must read a section of downsample in first!');
    disp('Use RSKreaddownsample, note that when dataset is small, downsample may not exist.')
    disp('If that is the case, try RSKreadthumbnail and RSKplotthumbnail.')
    handles = NaN;
    return
end



chanCol = [];
if ~strcmp(channel, 'all')
    channels = cellchannelnames(RSK, channel);
    for chan = channels
        chanCol = [chanCol getchannelindex(RSK, chan{1})];
    end
end

handles = channelsubplots(RSK, field, 'chanCol', chanCol);


end
