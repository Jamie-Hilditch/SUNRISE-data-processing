function handles = RSKplotthumbnail(RSK, varargin)

% RSKplotthumbnail - Plot summaries of logger data thumbnail.
%
% Syntax:  [handles] = RSKplotthumbnail(RSK, [OPTIONS])
% 
% Generates a summary plot of the thumbnail data in the RSK structure.
% ThumbnailData is usually a downsampled dataset of about 4000 points.
% Each time value has a maximum and a minimum data value so that all spikes
% are visible even though the dataset is down-sampled. 
% 
% Inputs:
%    [Required] - RSK - Structure containing the logger metadata and
%                       thumbnail.
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
%    RSKplotthumbmail(RSK);  
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



field = 'thumbnailData';
if ~isfield(RSK,field)
    disp('You must read a section of thumbnailData in first!');
    disp('Use RSKreadthumbnail...')
    disp('Note that thumbnail data may not be supported in future release.')
    disp('We recommend users to use RSKreaddownsample and RSKplotdownsample.')
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
