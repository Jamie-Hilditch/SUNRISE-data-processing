function data = rose_parse_deployment(dir_in, varargin)

    % dir_in = '/home/dw/Data/gdrive/ROSE/Deployments/202106_GoM/UBOX01_deploy_20210702_141010/';
    % progress = [];

    % Parse optional inputs
    p = inputParser;
    addParameter(p,'progress',[],@(x) isa(x,'matlab.ui.dialog.ProgressDialog'));
    parse(p,varargin{:});
    progress = p.Results.progress;
    parse_nuc_timestamps = true;

    % Find files in deployment folder
    files.adcp = dir(fullfile(dir_in,'ADCP','*ADCP_timestamped*.bin'));
    files.gps = dir(fullfile(dir_in,'GPS','GPS_*.log'));
    files.imu = dir(fullfile(dir_in,'IMU','IMU_timestamped*.bin'));
    isNortek = all(startsWith({files.adcp.name},'Nortek'));

    % Parse data
    data = struct('adcp',[],'gps',[],'imu',[],'stats',struct());

    fprintf('ADCP...')
    if ~isempty(files.adcp)
        if ~isempty(progress)
            if isNortek
                data.adcp = parse_nortek_adcp(files.adcp,parse_nuc_timestamps,...
                                              'progress',progress);
            else
                data.adcp = parse_adcp(files.adcp,parse_nuc_timestamps,...
                                       'progress',progress);
            end
        else
            if isNortek
                data.adcp = parse_nortek_adcp(files.adcp,parse_nuc_timestamps);
            else
                data.adcp = parse_adcp(files.adcp,parse_nuc_timestamps);
            end
        end
        if ~isempty(data.adcp)
            % data.adcp = fix_nuc_adcp_timestamps(data.adcp);
            data.stats.adcp_pings = length(data.adcp.time);
            data.stats.adcp_config = data.adcp.config;
            data.stats.adcp_fields = structfun(@size,data.adcp,'uni',false);
        end
    end
    fprintf(' Done\n')
    pause(1)

    fprintf('GPS...')
    if ~isempty(files.gps)
        if ~isempty(progress)
            data.gps = parse_gps(files.gps,'progress',progress);
        else
            data.gps = parse_gps(files.gps);
        end
        if ~isempty(data.gps)
            gps_fld_names = setdiff(fields(data.gps),'files');
            for i = 1:length(gps_fld_names)
                data.stats.gps_msg_count.(gps_fld_names{i}) = length(data.gps.(gps_fld_names{i}).lnum);
            end
        end
    end
    fprintf(' Done\n')
    pause(1)

    fprintf('IMU...')
    if ~isempty(files.imu)
        if ~isempty(progress)
            data.imu = parse_imu(files.imu,parse_nuc_timestamps,...
                                 'progress',progress);
        else
            data.imu = parse_imu(files.imu,parse_nuc_timestamps);
        end
        if ~isempty(data.imu)
            imuflds = setdiff(fields(data.imu),{'units'});
            for f = 1:length(imuflds)
                data.stats.imu_msg_count.(imuflds{f}) = ...
                    length(data.imu.(imuflds{f}).nuc_time);
                data.stats.imu_subflds.(imuflds{f}) = fields(data.imu.(imuflds{f}));
            end
        end
    end
    fprintf(' Done\n')
    pause(1)

    % Report total file sizes
    ftypes = fields(files);
    bytes_total = 0;
    for i = 1:length(ftypes)
        b = sum([files.(ftypes{i}).bytes]);
        data.stats.(sprintf('bytes_%s',ftypes{i})) = b;
        bytes_total = bytes_total + b;
    end
    data.stats.bytes_total = bytes_total;

end
