function RSK = RSKgetprofiles(RSK)

%RSKgetprofiles - Find the profiles start and end times.
%
% Syntax:  [RSK] = RSKgetprofiles(RSK)
% 
% Finds the profiles start and end times by first looking at the region
% table (Ruskin generated profiles) then at the events table (logger
% generated profiles). 
%
% Inputs: 
%    RSK - Structure, with profile events.
%
% Outputs:
%    RSK - Structure containing the logger metadata and thumbnails
%          including profile metadata.
%
% See also: RSKopen, RSKfindprofiles, RSKreadprofiles.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-22

if isfield(RSK, 'profiles')
    error('Profiles are already found, get data using RSKreadprofiles.m');
end



RSK = readregionprofiles(RSK);

if ~isfield(RSK, 'profiles')
    RSK = readeventsprofiles(RSK);
end

end



        
