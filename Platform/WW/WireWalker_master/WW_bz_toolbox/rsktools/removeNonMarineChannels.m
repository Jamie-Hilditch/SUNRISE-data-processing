function [RSK, isDerived] = removenonmarinechannels(RSK)

%REMOVENONMARINECHANNELS - Remove hidden or derived channels.
%
% Syntax:  [RSK, isDerived] = REMOVENONMARINECHANNELS(RSK)
%
% Removes the hidden or derived channels from the channels table and
% returns a logical index vector. They are also removed from
% instrumentChannels if the field exists. 
%
% Inputs:
%    RSK - Structure
%
% Outputs:
%    RSK - Structure with only marine channels.
%
%    isDerived - Logical index describing which channels are non-marine.
%
% See also: RSKopen, readheaderfull.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-21

p = inputParser;
addRequired(p, 'RSK', @isstruct);
parse(p, RSK)

RSK = p.Results.RSK;

if ~(strcmp(RSK.dbInfo(end).type, 'EPdesktop') || strcmp(RSK.dbInfo(end).type, 'skinny'))
    if iscompatibleversion(RSK, 1, 8, 9) && ~strcmp(RSK.dbInfo(end).type, 'EP')
        if logical(RSK.toolSettings.rhc)
            isDerived = logical([RSK.instrumentChannels.channelStatus] == 4);% derived channels have a '4' channelStatus
        else
            isDerived = logical([RSK.instrumentChannels.channelStatus]);% hidden and derived channels have a non-zero channelStatus
        end
        RSK.instrumentChannels(isDerived) = [];
    else
        results = doSelect(RSK, 'select isDerived from channels');
        isDerived = logical([results.isDerived])'; 
    end
else
    isDerived = false(length(RSK.channels));
end


if length(RSK.channels) == length(isDerived)
    RSK.channels(isDerived) = [];
end


end


