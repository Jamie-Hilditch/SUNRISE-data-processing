%% check_proc_methods.m
% Usage: methods = check_proc_methods(data)
% Description: Return a list of compatible processing methods based on present
%              data fields.
% Inputs: data - structure containing ADCP, GPS, and IMU substructures.
% Outputs: methods - structure containing names of valid processing methods.
%
% Calling this function with an empty data input and the 'all' parameter will
% return a list of all known methods. This can be used for populating lists or
% validating user input. For example:
%     proc_methods = check_proc_methods([],'all','adcp_orientation')
% returns a cell array of all adcp orientation methods.
function proc_methods = check_proc_methods(data,varargin)

% First list all available methods. This must be kept up-to-date as more methods
% are added.
proc_methods = struct();
proc_methods.vessel_vel = {...
    'GPRMC groundspeed and course';
    'GPRMC distance/time';
    'GPGGA distance/time';
    'GPZDA distance/time';
    'None'};
proc_methods.adcp_orientation = {...
    'HEHDT heading (GPRMC time) + offset';
    'HEHDT heading (GPGGA time) + offset';
    'HEHDT heading (GPZDA time) + offset';
    'Nortek AHRS + Hemisphere + offset';
    'Nortek AHRS + offset';
    'ADCP compass + offset';
    'ADCP compass';
    'GPRMC course + offset';
    'GPVTG course + offset'
    'Offset only'};

% If specified, return all methods (i.e. use this function to generate a list of
% all methods)
p = inputParser;
addParameter(p,'all',[],@(x) ismember(lower(x),fields(proc_methods)));
parse(p,varargin{:});
if ~isempty(p.Results.all)
    proc_methods = proc_methods.(p.Results.all);
    return
end

proc_method_flds = fields(proc_methods);
isvalid = struct();
for i = 1:length(proc_method_flds)
    isvalid.(proc_method_flds{i}) = false(size(proc_methods.(proc_method_flds{i})));
end

% Fill vessel velocity and adcp orientation methods based on GPS fields
if ~isfield(data,'gps') || isempty(data.gps)
    gps_flds = {};
else
    gps_flds = fields(data.gps);
end
for i = 1:length(gps_flds)
    switch upper(gps_flds{i})
      case 'GPRMC'
        isvalid.vessel_vel(find(strcmp('GPRMC distance/time', proc_methods.vessel_vel))) = true;
        isvalid.vessel_vel(find(strcmp('GPRMC groundspeed and course', proc_methods.vessel_vel))) = true;
        isvalid.adcp_orientation(find(strcmp('GPRMC course + offset', proc_methods.adcp_orientation))) = true;
      case {'GPGGA','GPZDA'}
        idx = find(strcmp(sprintf('%s distance/time',gps_flds{i}), proc_methods.vessel_vel));
        isvalid.vessel_vel(idx) = true;
        idx = find(strcmp(sprintf('%s course + offset',gps_flds{i}), proc_methods.adcp_orientation));
        isvalid.adcp_orientation(idx) = true;
      case 'HEHDT'
        for ii = 1:length(gps_flds)
            switch upper(gps_flds{ii})
              case {'GPRMC','GPGGA','GPZDA'}
                idx = find(strcmp(sprintf('%s heading (%s time) + offset',gps_flds{i},gps_flds{ii}),...
                                  proc_methods.adcp_orientation));
                 isvalid.adcp_orientation(idx) = true;
            end
        end
    end
end

% Need HEHDT & GPRMC for Nortek AHRS + Hemisphere method
if all(ismember({'HEHDT','GPRMC'},gps_flds))
    idx = find(strcmp('Nortek AHRS + Hemisphere + offset',proc_methods.adcp_orientation));
    isvalid.adcp_orientation(idx) = true;
    idx = find(strcmp('Nortek AHRS + offset',proc_methods.adcp_orientation));
    isvalid.adcp_orientation(idx) = true;
end

% The 'none' and 'offset only' method should always be valid
idx = find(strcmp('none',lower(proc_methods.vessel_vel)));
isvalid.vessel_vel(idx) = true;
idx = find(strcmp('offset only',lower(proc_methods.adcp_orientation)));
isvalid.adcp_orientation(idx) = true;

% Check that ADCP has heading field
if isfield(data,'adcp') && ~isempty(data.adcp)
    if isfield(data.adcp,'heading')
        idx = find(strcmp('ADCP compass', proc_methods.adcp_orientation));
        isvalid.adcp_orientation(idx) = true;
        idx = find(strcmp('ADCP compass + offset', proc_methods.adcp_orientation));
        isvalid.adcp_orientation(idx) = true;
    end
end

% Remove invalid fields from proc methods
for i = 1:length(proc_method_flds)
    proc_methods.(proc_method_flds{i}) = ...
        proc_methods.(proc_method_flds{i})(isvalid.(proc_method_flds{i}));
end

% TODO: Implement IMU-based methods
