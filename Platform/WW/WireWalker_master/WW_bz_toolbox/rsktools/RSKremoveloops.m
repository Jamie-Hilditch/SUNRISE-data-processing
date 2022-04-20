function [RSK, flagidx] = RSKremoveloops(RSK, varargin)

% RSKremoveloops - Remove data exceeding a threshold profiling rate and
% with reversed pressure (loops).
%
% Syntax:  [RSK, flagidx] = RSKremoveloops(RSK, [OPTIONS])
% 
% Identifies and flags data obtained when the logger vertical profiling
% speed falls below a threshold value or when the logger reversed the
% desired cast direction (forming a loop). The flagged data is replaced 
% with NaNs.  All logger channels except depth are affected.    
% 
% Differenciates depth to estimate the profiling speed. The depth channel
% is first smoothed with a 3-point running average to reduce noise. 
% 
% Inputs:
%   [Required] - RSK - RSK structure with logger data and metadata
%
%   [Optional] - profile - Profile number. Defaults to all profiles.
%
%                direction - 'up' for upcast, 'down' for downcast, or
%                      'both' for all. Defaults to all directions
%                       available.
% 
%                threshold - Minimum speed at which the profile must
%                       be taken. Defaults to 0.25 m/s.
%
% Outputs:
%    RSK - Structure with data filtered by threshold profiling speed and
%          removal of loops.
%
%    flagidx - Index of the samples that are filtered.
%
% Example: 
%    RSK = RSKopen('file.rsk');
%    RSK = RSKreadprofiles(RSK);
%    RSK = RSKremoveloops(RSK);
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-10-17

validDirections = {'down', 'up', 'both'};
checkDirection = @(x) any(validatestring(x,validDirections));

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addParameter(p, 'profile', [], @isnumeric);
addParameter(p, 'direction', [], checkDirection);
addParameter(p, 'threshold', 0.25, @isnumeric);
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
profile = p.Results.profile;
direction = p.Results.direction;
threshold = p.Results.threshold;



try
    Dcol = getchannelindex(RSK, 'Depth');
catch
    error('RSKremoveloops requires a depth channel to calculate velocity (m/s). Use RSKderivedepth...');
end



castidx = getdataindex(RSK, profile, direction);
for ndx = castidx
    d = RSK.data(ndx).values(:,Dcol);
    depth = runavg(d, 3, 'nan');
    time = RSK.data(ndx).tstamp;

    velocity = calculatevelocity(depth, time);
    
    if getcastdirection(depth, 'up')
      flag = (velocity > -threshold);
      cm = cummin(depth);
      flag((depth - cm) > 0) = true;
    else
      flag = velocity < threshold; 
      cm = cummax(depth);
      flag((depth - cm) < 0) = true;
    end
    
    flagChannels = ~strcmpi('Depth', {RSK.channels.longName});    
    RSK.data(ndx).values(flag,flagChannels) = NaN;
    flagidx(ndx).index = find(flag);
end



logdata = logentrydata(RSK, profile, direction);
logentry = ['Samples measured at a profiling velocity less than ' num2str(threshold) ' m/s were replaced with NaN on ' logdata '.'];

RSK = RSKappendtolog(RSK, logentry);

end
