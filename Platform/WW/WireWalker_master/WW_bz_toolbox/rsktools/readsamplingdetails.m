function RSK = readsamplingdetails(RSK)

%READSAMPLINGDETAILS - Read the sampling details of a file.
%
% Syntax:  [RSK] = READSAMPLINGDETAILS(RSK)
%
% Reads the table that contains the sampling detail, this depends on the
% mode of the file. Possible modes are ddsampling, directional, fetching or
% continuous. The change happened in RSK schema v1.13.8. 
%
% Inputs:
%    RSK - Structure containing some logger metadata.
%
% Output:
%    RSK - Structure containing previously present logger metadata as well
%          as sampling details.
%
% See also: readheaderfull, RSKopen.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-21

mode = RSK.schedules.mode;

if iscompatibleversion(RSK, 1, 13, 8)
    if strcmpi(mode, 'ddsampling')
        modetable = 'directional';
    elseif strcmpi(mode, 'fetching')
        return
    else 
        modetable = mode;
    end

    RSK.(modetable) = doSelect(RSK, ['select * from ' modetable]);
end

end