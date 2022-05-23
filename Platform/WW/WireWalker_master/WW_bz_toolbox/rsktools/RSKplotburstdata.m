function handles = RSKplotburstdata(RSK, varargin)

%RSKplotburstdata - Plot summaries of logger burst data.
%
% Syntax:  [handles] = RSKplotburstdata(RSK, [OPTIONS])
% 
% Generates a plot, similar to the thumbnail plot, only using the
% full 'burstData' that you read in, rather than just the thumbnail
% view.
% 
% Inputs:
%    [Required] - RSK - Structure containing the logger metadata and
%                       burstData.
%
%    [Optional] - channel - Longname of channel to plots, can be multiple
%                       in a cell, if no value is given it will plot all
%                       channels. 
%
% Output:
%     handles - Line object of the plot.
%
% Example: 
%    RSK = RSKopen('sample.rsk');  
%    RSK = RSKreadburstdata(RSK, 'channel', {'Conductivity', 'Temperature', 'Pressure'});  
%    RSKplotburstdata(RSK);  
%
% See also: RSKreadburstdata, RSKplotthumbnail, RSKplotdata.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-22

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addParameter(p, 'channel', 'all');
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
channel = p.Results.channel;



field = 'burstData';
if ~isfield(RSK, field)
    disp('You must read a section of burst data in first!');
    disp('Use RSKreadburstdata...')
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