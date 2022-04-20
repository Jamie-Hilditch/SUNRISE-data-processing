function RSK = readheaderlive(RSK)

%READHEADERLIVE - Read tables that are populated in a 'live' file.
%
% Syntax:  [RSK] = READHEADERLIVE(RSK)
%
% Opens the non-standard populated tables of RSK 'live' files. These tables
% are appSettings, instrumentsChannels, and parameters. If data is
% available, it will open the parameterKeys and thumbnailData tables.  
%
% Note: Only marine channels will be displayed.
%
% Inputs:
%    RSK - Structure of 'live' file opened using RSKopen.m.
%
% Outputs:
%    RSK - Structure containing logger metadata and thumbnail.
%
% See also: RSKopen.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-07-10

p = inputParser;
addRequired(p, 'RSK', @isstruct);
parse(p, RSK)

RSK = p.Results.RSK;

%% Tables that are definitely in 'live'
RSK.appSettings = doSelect(RSK, 'select * from appSettings');

RSK = readparameters(RSK);

RSK = readsamplingdetails(RSK);



%% Tables that may or may not be in 'live'
tables = doSelect(RSK, 'SELECT name FROM sqlite_master WHERE type="table"');

if any(strcmpi({tables.name}, 'geodata'))
    RSK = RSKreadgeodata(RSK);
end

if any(strcmpi({tables.name}, 'thumbnailData'))
    RSK = RSKreadthumbnail(RSK);
end

RSK = RSKreaddownsample(RSK);

end
