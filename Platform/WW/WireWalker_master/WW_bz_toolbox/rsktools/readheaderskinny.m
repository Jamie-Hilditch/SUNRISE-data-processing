function RSK = readheaderskinny(RSK)

%READHEADERSKINNY - Read populated tables of a 'skinny' file.
%
% Syntax:  [RSK] = READHEADERSKINNY(RSK)
%
% Opens the non-standard populated tables of 'skinny' files. Only to be
% used by RSKopen.m. If metadata is available, it will open geodata.
%
% Note: The data is stored in raw bin file, open this file in Ruskin first
% to read the data. 
%
% Inputs:
%    RSK - Structure of 'skinny' file opened using RSKopen.m.
%
% Outputs:
%    RSK - Structure containing logger metadata and thumbnail.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-21

tables = doSelect(RSK, 'SELECT name FROM sqlite_master WHERE type="table"');

if any(strcmpi({tables.name}, 'geodata'))
    RSK = RSKreadgeodata(RSK);
end

end