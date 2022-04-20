function [RSK] = RSKderivesalinity(RSK)

%RSKderivesalinity - Calculate practical salinity.
%
% Syntax:  [RSK] = RSKderivesalinty(RSK)
% 
% Derives salinity using the TEOS-10 GSW toolbox
% (http://www.teos-10.org/software.htm). The result is added to the
% RSK data structure, and the channel list is updated. If salinity is
% already in the RSK data structure (i.e., from Ruskin), it will be
% overwritten by RSKderivesalinity.
%
% Inputs: 
%    RSK - Structure containing the logger metadata and data.         
%
% Outputs:
%    RSK - Updated structure containing practical salinity.
%
% See also: RSKcalculateCTlag.
%
% Author: RBR Ltd. Ottawa ON, Canada
% email: support@rbr-global.com
% Website: www.rbr-global.com
% Last revision: 2017-07-02

if isempty(which('gsw_SP_from_C'))
    error('RSKtools requires TEOS-10 GSW toolbox to derive salinity. Download it here: http://www.teos-10.org/software.htm');
end
    


Ccol = getchannelindex(RSK, 'Conductivity');
Tcol = getchannelindex(RSK, 'Temperature');



RSK = addchannelmetadata(RSK, 'Salinity', 'PSU');
Scol = getchannelindex(RSK, 'Salinity');
[RSKsp, SPcol] = getseapressure(RSK);



castidx = getdataindex(RSK);
for ndx = castidx
    salinity = gsw_SP_from_C(RSK.data(ndx).values(:, Ccol), RSK.data(ndx).values(:, Tcol), RSKsp.data(ndx).values(:,SPcol));
    RSK.data(ndx).values(:,Scol) = salinity;
end

logentry = ('Practical Salinity derived using TEOS-10 GSW toolbox.');
RSK = RSKappendtolog(RSK, logentry);

end



