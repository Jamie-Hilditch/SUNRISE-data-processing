function data = proc_pressure_cal(data,cfg)

if isfield(cfg,'zero_pressure_interval') && ~isempty(cfg.zero_pressure_interval)
    % Sample onto pressure calibration interval
    pcalgrid = proc_grid_init(data,cfg,cfg.zero_pressure_interval);
    % Compute pressure offsets
    disp(sprintf('Calibrating pressure data over interval: %s,%s',...
                 datestr(cfg.zero_pressure_interval(1)),...
                 datestr(cfg.zero_pressure_interval(2))));
    for i = 1:length(pcalgrid.pos)
        if isfield(data{i},'P') && ~all(isnan(pcalgrid.P(i,:)));
            p0 = nanmean(pcalgrid.P(i,:));
            if ~isnan(p0);
                data{i}.P = data{i}.P - p0;
                disp(sprintf('Removed %.2fdbar pressure offset from %s',...
                             p0,data{i}.sn))
            end
        end
    end
end

