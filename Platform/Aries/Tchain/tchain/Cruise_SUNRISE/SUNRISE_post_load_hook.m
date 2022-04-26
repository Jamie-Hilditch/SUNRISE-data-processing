function [data, cfg] = SUNRISE_post_load_hook(data,cfg)

switch cfg.vessel
  case 'Aries'
    switch cfg.name
      case 'deploy_20210629'
        % This sensor had bad timestamps, weird pressure values... possibly
        % salvageable for 1st half of deployment.
        i = strcmp('60702',{cfg.sensors.sn});
        % dn_old = data{i}.dn;
        % data{i}.dn = data{i}.dn(1) + [0:length(data{i}.dn)-1]*1/16/86400;
        data = data(~i);
        cfg.sensors = cfg.sensors(~i);
    end

end
