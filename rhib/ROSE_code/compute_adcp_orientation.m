function [orientation, angles, vessel_heading] = compute_adcp_orientation(data,proc_method,varargin)

offset = 0;
if contains(lower(proc_method),'offset')
    p = inputParser;
    addRequired(p,'offset',@(x) isnumeric(x))
    parse(p,varargin{:});
    offset = p.Results.offset;
end
vessel_heading = 0;

switch proc_method
  case 'Nortek AHRS + offset'
    % Use the Nortek's internal AHRS to get data in a flat coordinate frame,
    % then rotate to ENU using Hemisphere heading.
    rot_order = 'ZYX';
    angles = [90-data.adcp.heading(:)+offset, -data.adcp.pitch(:), data.adcp.roll(:)];
    orientation = -quaternion(angles,'eulerd',rot_order,'frame');
  case 'Nortek AHRS + Hemisphere + offset'
    % Use the Nortek's internal AHRS to get data in a flat coordinate frame,
    % then rotate to ENU using Hemisphere heading.
    h = gps_line_interp(data.gps,'GPRMC','HEHDT','head','angular');
    [~,iu] = unique(data.gps.GPRMC.dn);
    h = gps_interp_heading(data.gps.GPRMC.dn(iu),h(iu),data.adcp.nuc_time);
    rot_order = 'ZYX';
    angles = [90-h(:)+offset, -data.adcp.pitch(:), data.adcp.roll(:)];
    orientation = -quaternion(angles,'eulerd',rot_order,'frame');
    vessel_heading = mod(h,360);
  case 'HEHDT heading (GPRMC time) + offset'
    % Use internal pitch/roll measurements, but use HEHDT heading and a mounting
    % offset for the final rotation about the vertical axis. Get HEHDT time from
    % nearby GPRMC timestamps.
    h = gps_line_interp(data.gps,'GPRMC','HEHDT','head','angular');
    [~,iu] = unique(data.gps.GPRMC.dn);
    h = gps_interp_heading(data.gps.GPRMC.dn(iu),h(iu),data.adcp.nuc_time);
    angles = [-data.adcp.roll(:), -data.adcp.pitch(:), h(:)+offset];
    rot_order = 'YXZ';
    orientation = quaternion(angles,'eulerd', rot_order, 'frame');
  case 'HEHDT heading (GPGGA time) + offset'
    % Use internal pitch/roll measurements, but use HEHDT heading and a mounting
    % offset for the final rotation about the vertical axis. Get HEHDT time from
    % nearby GPGGA timestamps.

    % Fix GPGGA timestamps for use as a time source
    dn0 = regexp(data.gps.files{1},'GPS_(\d+).log','tokens');
    dn0 = datenum(dn0{1}{1},'yyyymmdd');
    data.gps.GPGGA.dn = data.gps.GPGGA.dn - floor(data.gps.GPGGA.dn) + dn0;
    idx = find(diff(data.gps.GPGGA.dn) < 0) + 1;
    for ii = 1:length(idx)
        data.gps.GPGGA.dn(idx(ii):end) = data.gps.GPGGA.dn(idx(ii):end) + 1;
    end
    h = gps_line_interp(data.gps,'GPGGA','HEHDT','head','angular');
    [~,iu] = unique(data.gps.GPGGA.dn);
    h = gps_interp_heading(data.gps.GPGGA.dn(iu),h(iu),data.adcp.nuc_time);
    angles = [-data.adcp.roll(:), -data.adcp.pitch(:), h(:)+offset];
    rot_order = 'YXZ';
    orientation = quaternion(angles,'eulerd', rot_order, 'frame');
  case 'HEHDT heading (GPZDA time) + offset'
    % Use internal pitch/roll measurements, but use HEHDT heading and a mounting
    % offset for the final rotation about the vertical axis. Get HEHDT time from
    % nearby GPZDA timestamps.
    h = gps_line_interp(data.gps,'GPZDA','HEHDT','head','angular');
    [~,iu] = unique(data.gps.GPZDA.dn);
    h = gps_interp_heading(data.gps.GPZDA.dn(iu),h(iu),data.adcp.nuc_time);
    angles = [-data.adcp.roll(:), -data.adcp.pitch(:), h(:)+offset];
    rot_order = 'YXZ';
    orientation = quaternion(angles,'eulerd', rot_order, 'frame');
  case {'ADCP compass','ADCP compass + offset'}
    % The rotation to earth coordinates is a frame rotation about the ADCP's Y
    % (roll), X (pitch), then Z (vertical) axes.
    % The heading offset is the compass orientation of beam 3 when the vessel is
    % facing due North.
    angles = [-data.adcp.roll(:), -data.adcp.pitch(:), data.adcp.heading(:) + offset];
    rot_order = 'YXZ';
    orientation = quaternion(angles,'eulerd', rot_order, 'frame');
  case 'GPRMC course + offset'
  case 'GPVTG course + offset'
  case 'Offset only'
    % Rotate to flattened ship coordinates
    angles = [-data.adcp.roll(:), -data.adcp.pitch(:), 0*data.adcp.heading(:) + offset];
    rot_order = 'YXZ';
    orientation = quaternion(angles,'eulerd', rot_order, 'frame');
  case 'Internal AHRS + Hemisphere'
    h = gps_line_interp(data.gps,'GPRMC','HEHDT','head','angular');
    [~,iu] = unique(data.gps.GPRMC.dn);
    h = gps_interp_heading(data.gps.GPRMC.dn(iu),h(iu),data.adcp.nuc_time);
    angles = [h(:)+offset, -data.adcp.pitch(:), data.adcp.roll(:)];
    rot_order = 'ZXY';
    orientation = quaternion(angles,'eulerd', rot_order, 'frame');
  case 'Instrument coordinates'
    % Rotate to flattened ship coordinates
    angles = [0*data.adcp.roll(:), 0*data.adcp.pitch(:), 0*data.adcp.heading(:)];
    rot_order = 'YXZ';
    orientation = quaternion(angles,'eulerd', rot_order, 'frame');
end
