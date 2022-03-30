function txt = rose_data_report(data)
    txt = {sprintf('Total Size: %.2fMB\n',data.stats.bytes_total/2^20)};
    flds = fields(data);
    for i = 1:length(flds)
        switch flds{i}
          % case 'stats'
          case 'adcp'
            txt{end+1} = sprintf('* ADCP Data: %.2fMB\n',data.stats.bytes_adcp/2^20);
            if isempty(data.adcp)
                txt{end+1}= sprintf('  No data recorded\n');
            else
            if isfield(data.adcp,'nuc_time')
                txt{end+1}= sprintf('  - External clock:\n    %s\n    %s',...
                                    datestr(data.adcp.nuc_time(1)),...
                                    datestr(data.adcp.nuc_time(end)));
            end
            txt{end+1}= sprintf('  - Internal clock:\n    %s\n    %s',...
                                datestr(data.adcp.time(1)),...
                                datestr(data.adcp.time(end)));
                txt{end+1}= sprintf('  - Serial: %d\n',data.adcp.config.serial_number);
                txt{end+1}= sprintf('  - Frequency: %d kHz\n',data.adcp.config.frequency);
                txt{end+1}= sprintf('  - %d %0.2fm depth cells\n',...
                                    data.adcp.config.n_cells,data.adcp.config.depth_cell_length);
                txt{end+1}= sprintf('  - %d recorded pings\n',length(data.adcp.time));
            end
          case 'gps'
            txt{end+1} = sprintf('* GPS Data: %.2fMB\n',data.stats.bytes_gps/2^20);
            if isempty(data.gps)
                txt{end+1}= sprintf('  No data recorded\n');
            else
                if isfield(data.gps,'GPRMC')
                    txt{end+1}= sprintf('  - GPS clock:\n    %s\n    %s',...
                                        datestr(data.gps.GPRMC.dn(1)),...
                                        datestr(data.gps.GPRMC.dn(end)));
                end
                gpsflds = setdiff(fields(data.gps),{'files'});
                for f = 1:length(gpsflds)
                    txt{end+1}= sprintf('  - %s: %d records\n',gpsflds{f},length(data.gps.(gpsflds{f}).lnum));
                end
            end
          case 'imu'
            txt{end+1} = sprintf('* IMU Data: %.2fMB\n',data.stats.bytes_imu/2^20);
            if isempty(data.imu)
                txt{end+1} = '  No data recorded\n';
            else
                imuflds = setdiff(fields(data.imu),{'units'});
                for f = 1:length(imuflds)
                    txt{end+1}= sprintf('  - %s: %d records\n',...
                                        imuflds{f},data.stats.imu_msg_count.(imuflds{f}));
                    subflds = fields(data.imu.(imuflds{f}));
                    for sf = 1:length(subflds)
                        txt{end+1}= sprintf('    - %s\n',subflds{sf});
                    end
                end
            end
        end
    end
end
