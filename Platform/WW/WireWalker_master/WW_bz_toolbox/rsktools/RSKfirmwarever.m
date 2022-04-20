function [v, vsnMajor, vsnMinor]  = RSKfirmwarever(RSK)

%RSKfirmwarever - Return the firmware version of the RSK file.
%
% Syntax:  [v, vsnMajor, vsnMinor] = RSKfirmwarever(RSK)
%
% Returns the most recent version of the firmware; the information is
% retrieved from 'instruments' fields for files older than v1.12.2 or
% 'deployments' for more recent files.
%
% Inputs:
%    RSK - Structure containing the logger metadata and thumbnail
%          returned by RSKopen.
%
% Output:
%    v - Lastest version of the firmware
%    vsnMajor - Latest version number of category major
%    vsnMinor - Latest version number of category minor
%    vsnPatch - Latest version number of category patch.
%
% See also: RSKver.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-22

if iscompatibleversion(RSK, 1, 12, 2)
    v = RSK.instruments.firmwareVersion;
else
    v = RSK.deployments.firmwareVersion;
end



vsn = textscan(v,'%s','delimiter','.');
vsnMajor = str2double(vsn{1}{1});
vsnMinor = str2double(vsn{1}{2});

end