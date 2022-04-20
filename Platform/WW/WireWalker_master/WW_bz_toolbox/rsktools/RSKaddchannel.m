function [RSK] = RSKaddchannel(RSK, newChan, channelName, units)

% RSKaddchannel - Add a new channel with defined channel name and
% units. If the new channel already exists in the RSK structure, it
% will overwrite the old one.
%
% Syntax:  [RSK] = RSKaddchannel(RSK, newChan, channelName, units)
% 
% Inputs: 
%    RSK - Structure containing the logger metadata and data. 
%    
%    newChan - Structure containing the data to be added.  The data for
%            the new channel must be stored in a field of newChan
%            called "values" (i.e., newChan.values).  If the data is
%            arranged as profiles in the RSK structure, then newChan
%            must be a 1xN array of structures of where N =
%            length(RSK.data).
% 
%    channelName - name of the added channel
%
%    units - unit of the added channel
%
% Outputs:
%    RSK - Updated structure containing the new channel.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2018-01-24

p = inputParser;
addRequired(p, 'RSK', @isstruct);
addRequired(p, 'newChan',@isstruct);
addRequired(p, 'channelName', @ischar);
addRequired(p, 'units', @ischar);
parse(p, RSK, newChan, channelName, units)

RSK = p.Results.RSK;
newChan = p.Results.newChan;
channelName = p.Results.channelName;
units = p.Results.units;

RSK = addchannelmetadata(RSK, channelName, units);
Ncol = getchannelindex(RSK, channelName);
castidx = getdataindex(RSK);
    
for ndx = castidx
   if ~isequal(size(newChan(ndx).values), size(RSK.data(ndx).tstamp));
       error('The dimensions of the new channel data structure must be consistent with RSK structure.')
   else
       RSK.data(ndx).values(:,Ncol) = newChan(ndx).values(:);
   end
end

logentry = [channelName ' (' units ') added to data table by RSKaddchannel'];
RSK = RSKappendtolog(RSK, logentry);


end