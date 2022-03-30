function adcp = rose_proc_adcp(data_raw,adcp_orientation_method,vessel_velocity_method)

%% Parse input arguments
p = inputParser;
addRequired(p,'adcp_orientation_method',@check_adcp_orientation_method)
addRequired(p,'vessel_velocity_method',@check_vessel_vel_method)
parse(p,adcp_orientation_method,vessel_velocity_method);
adcp_orientation_method = p.Results.adcp_orientation_method;
vessel_velocity_method = p.Results.vessel_velocity_method;

adcp = data_raw.adcp;

% FIXME correct nuc timestamps by interpolating using only timestamps in the
% correct decade... Probably a better way to do this.
idx = adcp.nuc_time > datenum([2020 0 0 0 0 0]);
adcp.nuc_time = interp1(find(idx),adcp.nuc_time(idx),1:length(adcp.nuc_time));

%% Compute ADCP orientation and vessel velocity
if isstr(adcp_orientation_method)
    [adcp_orientation, angles, vessel_heading] = ...
        compute_adcp_orientation(data_raw,adcp_orientation_method);
else
    [adcp_orientation, angles, vessel_heading] = ...
        compute_adcp_orientation(data_raw,adcp_orientation_method{:});
    adcp.computed_heading = mod(angles(:,3),360);
end

if isstr(vessel_velocity_method)
    vessel_vel0 = compute_vessel_vel(data_raw,vessel_velocity_method);
else
    vessel_vel0 = compute_vessel_vel(data_raw,vessel_velocity_method{:});
end

adcp.vel = vel_beam2earth(adcp.vel,adcp.config,adcp_orientation);

% Interpolate vessel velocity to ADCP time,
% Remove vessel motion from East/West and North/South components of ADCP signal
vessel_vel = nan(length(adcp.time),2);

%% Filtering

% Remove data below 80% of bottom track depth
[bd2 cd2] = meshgrid(0.80*cosd(25)*nanmean(adcp.bt_range), adcp.cell_depth');
mask = +permute(cd2 <= bd2, [1 3 2]);
mask(mask==0) = nan;
adcp.vel = adcp.vel .* mask;

% Remove data where error vel exceeds some threshold
% evel_max = 0.1;
% mask = +(adcp.vel(:,4,:) < evel_max);
% mask(mask==0) = nan;
% adcp.vel(:,1:3,:) = adcp.vel(:,1:3,:) .* mask;

% Interpolate vessel vel to ADCP clock
for i = 1:2
    vessel_vel(:,i) = interp1(vessel_vel0.time, vessel_vel0.vel(:,i), adcp.nuc_time(:));
    % Don't compensate for vessel motion if BT has already been removed
    if ~(isfield(adcp.processing,'bt_removed_from_vel') && adcp.processing.bt_removed_from_vel)
        adcp.vel(:,i,:) = squeeze(adcp.vel(:,i,:)) + repmat([vessel_vel(:,i)]', adcp.config.n_cells,1);
    end
end

adcp.processing.adcp_orientation_method = adcp_orientation_method;
adcp.processing.vessel_velocity_method = vessel_velocity_method;
adcp.vessel_u = vessel_vel(:,1);
adcp.vessel_v = vessel_vel(:,2);
if ~isempty(vessel_heading)
    adcp.vessel_heading = vessel_heading;
end
adcp.vessel_lat = interp1(vessel_vel0.time, vessel_vel0.lat, adcp.nuc_time);
adcp.vessel_lon = interp1(vessel_vel0.time, vessel_vel0.lon, adcp.nuc_time);
