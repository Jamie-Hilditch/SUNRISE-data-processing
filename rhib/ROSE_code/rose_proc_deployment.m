function proc = rose_proc_deployment(data_raw,adcp_orientation_method,vessel_velocity_method)

%% Parse input arguments
p = inputParser;
addRequired(p,'adcp_orientation_method',@check_adcp_orientation_method)
addRequired(p,'vessel_velocity_method',@check_vessel_velocity_method)
parse(p,adcp_orientation_method,vessel_velocity_method);
adcp_orientation_method = p.Results.adcp_orientation_method;
vessel_velocity_method = p.Results.vessel_velocity_method;

adcp_orientation = compute_adcp_orientation_method(data_raw,adcp_orientation_method{:});
vessel_vel = compute
