function samplingperiod = RSKsamplingperiod(RSK)

%RSKsamplingperiod - Returns the sampling period information.
%
% Syntax:  [samplingperiod] = RSKsamplingperiod(RSK)
%
% Returns the sampling period of the file.
%
% Inputs:
%    RSK - Structure containing the logger metadata and thumbnail.
%
% Output:
%    samplingperiod - In seconds.
%
% See also: RSKfirmwarever, RSKver.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-22

mode = RSK.schedules.mode;
if iscompatibleversion(RSK, 1, 13, 8)
    if strcmpi(mode, 'ddsampling')
        samplingperiod.fastThreshold = RSK.directional.fastThreshold/1000;
        samplingperiod.slowThreshold = RSK.directional.slowThreshold/1000;
    elseif strcmpi(mode, 'fetching')
        error('"Fetching" files do not have a sampling period');
    else 
        samplingperiod = RSK.(mode).samplingPeriod/1000;
    end
else
    samplingperiod = RSK.schedules.samplingPeriod/1000;
end

end