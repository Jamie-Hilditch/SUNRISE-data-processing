function RSK = readheaderEPdesktop(RSK)

%READHEADEREPDESKTOP - Read tables that are populated in a 'EPdesktop' file.
%
% Syntax:  [RSK] = READHEADEREPDESKTOP(RSK)
%
% Opens the non-standard populated tables of 'EPdesktop' files, this table 
% is thumbnailData. If metadata is present, it will also read the
% appSettings, parameters, parameterKeys and geodata tables.  
%
% Inputs:
%    RSK - Structure of 'EPdesktop' file opened using RSKopen.m.
%
% Outputs:
%    RSK - Structure containing logger metadata and thumbnail.
%
% See also: RSKopen.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-21

%% Tables that are definitely in 'EPdesktop'
RSK = RSKreadthumbnail(RSK);

RSK = readsamplingdetails(RSK);

if ~strcmpi(RSK.dbInfo(1).type, 'EasyParse')
    if ~strcmpi(RSK.dbInfo(1).type, 'skinny')
        RSK = readparameters(RSK);
    end
end



%% Tables that may or may not be in file
tables = doSelect(RSK, 'SELECT name FROM sqlite_master WHERE type="table"');

if any(strcmpi({tables.name}, 'geodata'))
    RSK = RSKreadgeodata(RSK);
end

if any(strcmpi({tables.name}, 'appSettings'))
    RSK.appSettings = doSelect(RSK, 'select * from appSettings');  
end

end

