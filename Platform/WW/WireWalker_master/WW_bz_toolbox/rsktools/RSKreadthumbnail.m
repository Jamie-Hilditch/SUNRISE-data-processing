function RSK = RSKreadthumbnail(RSK)

%RSKreadthumbnail - Read thumbnail data from an opened RSK file. 
%
% Syntax:  [RSK] = RSKreadthumbnail(RSK)
% 
% Reads thumbnail data from an opened RSK SQLite file, called from
% within RSKopen.
%
% Inputs:
%    RSK - Structure containing the logger metadata and thumbnails
%          returned by RSKopen.
%
% Output:
%    RSK - Structure containing previously present logger metadata as well
%          as thumbnailData.
%
% See also: RSKopen, RSKplotthumbnail.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-22

p = inputParser;
addRequired(p, 'RSK', @isstruct);
parse(p, RSK);

RSK = p.Results.RSK;

sql = 'select tstamp/1.0 as tstamp, * from thumbnailData order by tstamp';
results = doSelect(RSK, sql);
if isempty(results)
    return
end



results = removeUnusedDataColumns(results);
results = arrangedata(results);

results.tstamp = RSKtime2datenum(results.tstamp');

if ~strcmpi(RSK.dbInfo(end).type, 'EPdesktop')
    [~, isDerived] = removeNonMarineChannels(RSK);
    results.values = results.values(:,~isDerived);
end



RSK.thumbnailData = results;

end
