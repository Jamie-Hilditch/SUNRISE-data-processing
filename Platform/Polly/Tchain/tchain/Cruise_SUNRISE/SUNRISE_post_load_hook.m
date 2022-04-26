function [data, cfg] = SUNRISE_post_load_hook(data,cfg)

switch cfg.vessel

  case 'Polly'
    switch cfg.name
      case 'deploy_20210625'
        % This deployment had 1 sensor clock way out of sync. Clock on recovery
        % showed 20000102_162841 while computer UTC clock showed
        % 20210627_134659.
        idx = find(strcmp({cfg.sensors.sn},'207057'));
        bad_time = datenum([2000 01 02 16 28 41]);
        good_time = datenum([2021 06 27 13 46 59]);
        data{idx}.dn = data{idx}.dn + (good_time-bad_time);

        % RBR Solo 101179 was having issues and didn't record any data
        rm = strcmp({cfg.sensors.sn},'101179');
        data = data(~rm);
        cfg.sensors = cfg.sensors(~rm);

        % Fix start/end times temporarily
        % cfg.dn_range = [min(cellfun(@(c) c.dn(1),data)),...
        %                 max(cellfun(@(c) c.dn(end),data))];
    end


end
