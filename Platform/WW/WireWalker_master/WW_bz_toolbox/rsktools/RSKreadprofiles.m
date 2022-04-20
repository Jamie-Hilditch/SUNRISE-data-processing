function RSK = RSKreadprofiles(RSK, varargin)

% RSKreadprofiles - Read individual casts from RSK SQLite database or
% existing RSK.data field.
%
% Syntax:  [RSK] = RSKreadprofiles(RSK, [OPTIONS])
% 
% Reads profile, including upcasts, downcasts, or both from the events 
% contained in a .rsk file. Each cast is an element in the data field 
% matrix. The cast direction is indicated as 'up' or 'down' in 
% RSK.data.direction.
%
% The profile events are parsed from the events table using the
% following types (see RSKconstants.m):
%   33 - Begin upcast
%   34 - Begin downcast
%   35 - End of profile cast
%
% Inputs: 
%    [Required] - RSK - Structure containing the logger data read
%                       from the RSK file.
%
%    [Optional] - profile - Vector identifying the profile numbers to
%                       read. Can be used to read only a subset of all
%                       the profiles. Default is to read all the profiles. 
% 
%                 direction - 'up' for upcast, 'down' for downcast, or
%                       `both` for all. Default is 'both'.
%
% Outputs:
%    RSK - RSK structure containing individual casts as each element in the
%          data field.
%
% Examples:
%    rsk = RSKopen('profiles.rsk');
%
%    % read all profiles
%    rsk = RSKreadprofiles(rsk);
%    -OR-
%    % read selective upcasts
%    rsk = RSKreadprofiles(rsk, 'profile', [1 3 10], 'direction', 'up');
%
% See also: RSKfindprofiles, RSKplotprofiles.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-11-22

validDirections = {'down', 'up', 'both'};
checkDirection = @(x) any(validatestring(x,validDirections));

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addParameter(p, 'profile', [], @isnumeric);
addParameter(p, 'direction', 'both', checkDirection);
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
profile = p.Results.profile;
direction = {p.Results.direction};



if ~isfield(RSK, 'profiles') 
    error('No profiles in this RSK, try RSKreaddata or RSKfindprofiles');
end
if strcmpi(direction{1}, 'both')
    direction = {'down', 'up'};
end



alltstart = [];
alltend = [];
for dir = direction
    castdir = [dir{1} 'cast'];
    alltstart = [alltstart; RSK.profiles.(castdir).tstart];
    alltend = [alltend; RSK.profiles.(castdir).tend];
end
alltstart = sort(alltstart);
alltend = sort(alltend);

RSK.profiles.order = direction;
profilecast = size(RSK.profiles.order, 2);
if profilecast == 2 && (alltstart(1) == RSK.profiles.upcast.tstart(1))
    RSK.profiles.order = {'up', 'down'};
end



if ~isempty(profile)
    if max(profile) > length(alltstart)/profilecast
        disp('The profile selected is greater than the total amount of profiles in this file.');
        return
    end
    if profilecast == 2
        castidx = [(profile*2)-1 profile*2];
        castidx = sort(castidx);
    else
        castidx = profile;
    end
else
    castidx = 1:length(alltstart);
end
RSK.profiles.originalindex = castidx;

dir2fill = cell(length(castidx),1); % append data.direction to each cast
if size(RSK.profiles.order, 2) == 1
    dir2fill(:) = direction;
    pronum2fill = castidx;
else
    dir2fill(1:2:end) = RSK.profiles.order(1);
    dir2fill(2:2:end) = RSK.profiles.order(2);
    pronum2fill = reshape(repmat(castidx(1:length(castidx)/2), 2, 1),length(castidx),1);
end

k = 1;
data(length(castidx)).tstamp = [];
data(length(castidx)).values = [];
data(length(castidx)).direction = [];
data(length(castidx)).profilenumber = [];

for ndx = castidx
    
    if isfield(RSK, 'data')
        ind_start = (find(RSK.data.tstamp == alltstart(ndx)));
        ind_end = (find(RSK.data.tstamp == alltend(ndx)));
        
        if isempty(ind_start) || isempty(ind_end)
            tmp = RSKreaddata(RSK, 't1', alltstart(ndx), 't2', alltend(ndx));
            data(k).tstamp = tmp.data.tstamp;
            data(k).values = tmp.data.values;
        else
            data(k).tstamp = RSK.data.tstamp(ind_start:ind_end);
            data(k).values = RSK.data.values(ind_start:ind_end,:);
        end
        
    else
        tmp = RSKreaddata(RSK, 't1', alltstart(ndx), 't2', alltend(ndx));
        data(k).tstamp = tmp.data.tstamp;
        data(k).values = tmp.data.values;
    end

    data(k).direction = dir2fill{k};
    data(k).profilenumber = pronum2fill(k);
    k = k + 1;
    
end

if ~isfield(RSK, 'data'), RSK = readchannels(RSK); end
RSK.data = data;

end