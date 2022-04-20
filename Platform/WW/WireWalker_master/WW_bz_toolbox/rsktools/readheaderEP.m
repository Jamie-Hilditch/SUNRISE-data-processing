function RSK = readheaderEP(RSK)

%READHEADEREP - Read tables that are populated in a 'EasyParse' file.
%
% Syntax:  [RSK] = READHEADEREP(RSK)
%
% For the tables that could be in an 'EasyParse' file, checks if they are
% present and populated and adds them to the RSK structure. 
%
% Inputs:
%    RSK - Structure of 'EasyParse' file opened using RSKopen.m.
%
% Outputs:
%    RSK - Structure containing the logger metadata.
%
% See also: RSKopen.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-07-10

%% Tables that could be populated in 'EasyParse'
tables = doSelect(RSK, 'SELECT name FROM sqlite_master WHERE type="table"');

if any(strcmpi({tables.name}, 'geodata'))
    RSK = RSKreadgeodata(RSK);
end

end

