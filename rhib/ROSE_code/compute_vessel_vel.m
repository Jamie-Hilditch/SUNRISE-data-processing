% TODO: Add function description
% TODO: Add rotational velocity components
function vessel_vel = compute_vessel_vel(data,proc_method)

vessel_vel = struct();

switch proc_method
  case 'GPRMC groundspeed and course'
    % pick unique timestamps to avoid interpolation issues later
    [~,iu] = unique(data.gps.GPRMC.dn);
    vessel_vel.time = data.gps.GPRMC.dn(iu); % GPRMC timestamps
    % Compute EW and NS velocities
    vessel_vel.vel = zeros(length(iu),3);
    vessel_vel.vel(:,1) = data.gps.GPRMC.speed(iu) .* cosd(90-data.gps.GPRMC.course(iu)); % E/W vel
    vessel_vel.vel(:,2) = data.gps.GPRMC.speed(iu) .* sind(90-data.gps.GPRMC.course(iu)); % N/S vel
    vessel_vel.lat = data.gps.GPRMC.lat(iu);
    vessel_vel.lon = data.gps.GPRMC.lon(iu);

  case 'GPRMC distance/time'
    [vessel_vel.time,iu] = unique(data.gps.GPRMC.dn);
    vessel_vel.vel = zeros(length(iu),3);
    [vessel_vel(:,1), vessel_vel(:,2)] = gps_ltln2vel(...
        data.gps.GPRMC.lat(iu),...
        data.gps.GPRMC.lon(iu),...
        data.gps.GPRMC.dn(iu));
    vessel_vel.lat = data.gps.GPRMC.lat(iu);
    vessel_vel.lon = data.gps.GPRMC.lon(iu);

  case 'GPGGA distance/time'
    % Fix GPGGA timestamps for use as a time source
    dn0 = regexp(data.gps.files{1},'GPS_(\d+).log','tokens');
    dn0 = datenum(dn0{1}{1},'yyyymmdd');
    data.gps.GPGGA.dn = data.gps.GPGGA.dn + dn0;
    idx = find(diff(data.gps.GPGGA.dn) < 0) + 1;
    for ii = 1:length(idx)
        data.gps.GPGGA.dn(idx(ii):end) = data.gps.GPGGA.dn(idx(ii):end) + 1;
    end
    [vessel_vel.time,iu] = unique(data.gps.GPGGA.dn);
    vessel_vel.vel = zeros(length(iu),3);
    [vessel_vel.vel(:,1), vessel_vel.vel(:,2)] = gps_ltln2vel(...
        data.gps.GPGGA.lat(iu),...
        data.gps.GPGGA.lon(iu),...
        data.gps.GPGGA.dn(iu));
    vessel_vel.vel(:,3) = gradient(data.gps.GPGGA.alt)./(gradient(data.gps.GPGGA.dn)*86400);
    vessel_vel.lat = data.gps.GPGGA.lat(iu);
    vessel_vel.lon = data.gps.GPGGA.lon(iu);

  case 'GPZDA distance/time'
  case 'IMU GNSS'
    vesse_vel = struct();
    [vessel_vel.time,iu] = unique(data.imu.GNSS.nuc_time);
    vessel_vel.vel = zeros(length(iu),3);
    [vessel_vel.vel(:,1), vessel_vel.vel(:,2)] = gps_ltln2vel(...
        data.imu.GNSS.llh_position.latitude(iu),...
        data.imu.GNSS.llh_position.longitude(iu),...
        data.imu.GNSS.nuc_time(iu)');
    vessel_vel.vel(:,3) = gradient(data.imu.GNSS.nuc_time(iu))./...
        (gradient(data.imu.GNSS.llh_position.height_above_msl(iu)')*86400);
  case 'None'
    vessel_vel = struct();
    if isfield(data.adcp,'nuc_time')
        vessel_vel.time = data.adcp.nuc_time;
    else
        vessel_vel.time = data.adcp.time;
    end
    vessel_vel.vel = zeros(length(data.adcp.time),3);
    vessel_vel.lat = 0*data.adcp.time;
    vessel_vel.lon = 0*data.adcp.time;
end
