function castidx = getdataindex(RSK, varargin)

%GETDATAINDEX - Return the index of data elements requested.
%
% Syntax:  [castIdx] = GETDATAINDEX(RSK, [OPTIONS])
% 
% Selects the data elements that fulfill the requirements described by the 
% profile number and direction arguments.
%
% Inputs:
%   [Required] - RSK - Structure containing the logger data
%
%   [Optional] - profile - Profile number. Default is to use all profiles
%                      available.
% 
%                direction - Cast direction. Default is to use all
%                      directions available. 
%            
% Outputs:
%    castidx - Array containing the index of data's elements.
%
% See also: RSKplotprofile, RSKsmooth, RSKdespike.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-20

validationFcn = @(x) ischar(x) || isempty(x);

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addOptional(p, 'profile', [], @isnumeric);
addOptional(p, 'direction', [], validationFcn);
parse(p, RSK, varargin{:})

RSK = p.Results.RSK;
profile = p.Results.profile;
direction = p.Results.direction;


if size(RSK.data,2) == 1
    castidx = 1;
    if ~isempty(profile) && profile ~= 1  
        error('The profile requested is greater than the total amount of profiles in this RSK structure.');
    end 
    return
end

profilecast = size(RSK.profiles.order,2);
ndata = length(RSK.data);

if ~isempty(direction) && profilecast == 1 && ~strcmp(RSK.profiles.order, direction)
    error(['There is no ' direction 'cast in this RSK structure.']);
end



if isempty(profile) && isempty(direction)
    castidx = 1:ndata;
elseif ~isempty(profile)
    if max(profile) > ndata/profilecast
        error('The profile requested is greater than the total amount of profiles in this RSK structure.');
    end
    
    if profilecast == 2
        if isempty(direction) || strcmp(direction, 'both')
            castidx = [(profile*2)-1 profile*2];
            castidx = sort(castidx);
        elseif strcmp(RSK.profiles.order{1}, direction)
            castidx = (profile*2)-1;
        else
            castidx = profile*2;
        end
    else
        castidx = profile;
    end
else
    if profilecast == 2
        if strcmp(direction, 'both')
            castidx = 1:ndata;
        elseif strcmp(RSK.profiles.order{1}, direction)
            castidx = 1:2:ndata;
        else
            castidx = 2:2:ndata;
        end
    else
        castidx = 1:ndata;
    end
end
     
end

