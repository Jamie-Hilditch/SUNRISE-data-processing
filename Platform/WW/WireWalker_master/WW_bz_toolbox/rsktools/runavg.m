function out = runavg(in, windowLength, edgepad)

%RUNAVG - Smooth a time series using a boxcar filter.
%
% Syntax:  [out] = RUNAVG(in, windowLength, edgepad)
% 
% Performs a running average, also known as boxcar filter, of length
% windowLength over the time series. 
%
% Inputs:
%    in - Time series
%
%    windowLength - Length of the running median. It must be odd.
%
%    edgepad - Describes how the filter will act at the edges. Options
%         are 'mirror', 'zeroorderhold' and 'nan'. Default is 'mirror'.
%
% Outputs:
%    out - Smoothed time series.
%
% See also: RSKsmooth, RSKcalculateCTlag.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-21

if nargin == 2
    edgepad = 'mirror';
end

if mod(windowLength, 2) == 0
    error('windowLength must be odd');
end



padsize = (windowLength-1)/2;
inpadded = padseries(in, padsize, edgepad);



n = length(in);
out = NaN*in;
for ndx = 1:n
    out(ndx) = nanmean(inpadded(ndx:ndx+(windowLength-1)));
end

end