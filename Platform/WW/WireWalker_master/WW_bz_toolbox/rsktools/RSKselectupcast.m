function [RSK, upidx] = RSKselectupcast(RSK)

%RSKselectupcast - Select the data elements with decreasing pressure.
%
% Syntax:  [RSKup, upidx] = RSKselectupcast(RSK)
%
% Keeps only the upcasts in the RSK and returns the index of the upcasts from
% the input RSK structure. 
%
% Inputs:
%    RSK - Structure containing logger data.
%
% Outputs:
%    RSK - Structure only containing upcasts.
%
%    upidx - Index of the data fields from the input RSK structure that 
%          are upcasts.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-06-23

Pcol = getchannelindex(RSK, 'Pressure');
ndata = length(RSK.data);



upidx = NaN(1, ndata);
for ndx = 1:ndata
    pressure = RSK.data(ndx).values(:, Pcol);
    upidx(1, ndx) = getcastdirection(pressure, 'up');
end



if ~any(upidx ==1)
    disp('No upcasts in this RSK structure.');
    return;
end



RSK.profiles.originalindex = RSK.profiles.originalindex(logical(upidx));
RSK.profiles.order = {'up'};
RSK.data = RSK.data(logical(upidx));

end