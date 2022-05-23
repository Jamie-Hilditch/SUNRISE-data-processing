function handles = channelsubplots(RSK, field, varargin)

%CHANNELSUBPLOTS - Plot each channel specified in a different subplot.
%
% Syntax:  [handles] = CHANNELSUBPLOT(RSK, field, [OPTIONS])
% 
% Generates subplots and plots each channel in the chosen data element.
% If data has many fields and none are specified, the first one is
% selected.  
%
% Inputs:
%   [Required] - RSK - Structure created by a .rsk file
%
%                field - Source of the data to plot. Can be
%                      'burstData', 'thumbnailData', or 'data'.
%
%   [Optional] - chanCol - Column number of the channels to plot.
%                      Default is to plot all channels.
%
%                castidx - Data element that will be used to make
%                      the plot. The default is 1. Note: To compare data
%                      elements use RSKplotprofiles. 
%
% Outputs:
%    handles - Line object of the plot.
%
% See also: RSKplotthumbnail, RSKplotburstdata, RSKplotdata.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-19

validFields = {'burstData', 'thumbnailData', 'data','downsample'};
checkField = @(x) any(validatestring(x,validFields));

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addRequired(p, 'field', checkField);
addParameter(p, 'chanCol', [], @isnumeric);
addParameter(p, 'castidx', 1);
parse(p, RSK, field, varargin{:});

RSK = p.Results.RSK;
field = p.Results.field;
chanCol = p.Results.chanCol;
castidx = p.Results.castidx;



if isempty(chanCol)
    chanCol = 1:size(RSK.(field)(castidx).values,2);
end
numchannels = length(chanCol);

n = 1;
for chan = chanCol
    subplot(numchannels,1,n)
    handles(n) = plot(RSK.(field)(castidx).tstamp, RSK.(field)(castidx).values(:,chan),'-');
    title(RSK.channels(chan).longName);
    ylabel(RSK.channels(chan).units);
    ax(n)=gca;
    datetick('x');
    n = n+1;
end

linkaxes(ax,'x');
shg

end