function RSK = readregionprofiles(RSK)

%READREGIONPROFILES - Read profiles start and end times from regions table.
%
% Syntax:  [RSK] = READREGIONPROFILES(RSK)
%
% Reads in profiles start and end time by combining information in the
% region and regionCast tables and adds it to the RSK structure.
%
% Inputs:
%    RSK - Structure containing logger metadata.
%
% Outputs:
%    RSK - Structure containing populated profiles, if available.
%
% See also: RSKgetprofiles.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-21

tables = doSelect(RSK, 'SELECT name FROM sqlite_master WHERE type="table"');

if any(strcmpi({tables.name}, 'regionCast')) && any(strcmpi({tables.name}, 'region'))
    RSK.region = doSelect(RSK, 'select regionID, type, tstamp1/1.0 as tstamp1, tstamp2/1.0 as tstamp2 from region');
    RSK.regionCast = doSelect(RSK, 'select * from regionCast');
else
    return
end



if ~isempty(RSK.regionCast)
    if strcmpi(RSK.regionCast(1).type, 'down')
        firstdir = 'downcast';
        lastdir = 'upcast';
    else
        firstdir = 'upcast';
        lastdir = 'downcast';
    end
else
    RSK = rmfield(RSK, 'regionCast');
    RSK = rmfield(RSK, 'region');
    return
end



for ndx = 1:length(RSK.regionCast)/2
    nregionCast = (ndx*2)-1;
    
    regionID = RSK.regionCast(nregionCast).regionID;
    RSK.profiles.(firstdir).tstart(ndx,1) = RSKtime2datenum(RSK.region(regionID).tstamp1); 
    RSK.profiles.(firstdir).tend(ndx,1) = RSKtime2datenum(RSK.region(regionID).tstamp2);
    
    regionID2 = RSK.regionCast(nregionCast+1).regionID;
    RSK.profiles.(lastdir).tstart(ndx,1) = RSKtime2datenum(RSK.region(regionID2).tstamp1);
    RSK.profiles.(lastdir).tend(ndx,1) = RSKtime2datenum(RSK.region(regionID2).tstamp2);
end
end
