function [data, cfg] = SUNRISE_post_load_hook(data,cfg)

% Note: I don't recommend modifying this in the version we're using as a group.
% The purpose of these hook functions is to arbitrarily modify the data at any
% stage in processing to do custom, user-specific tasks.
idx = find(strcmp({cfg.sensors.sn},'60523'));
if ~isempty(idx)
    tmp = data{idx}.p;
    data{idx}.p = data{idx}.t;
    data{idx}.t = tmp;
    
    data{idx}.s = gsw_SP_from_C(data{idx}.c, data{idx}.t, 0);
end

% idx2 = find(strcmp({cfg.name}, 'deploy_20210623'));
% if ~isempty(idx2)
%     idx3 = find(strcmp({cfg.sensors.sn},'60559'));
%     if ~isempty(idx)
%         data{idx3}.dn(1766064:end) = data{idx3}.dn(1766064:end) - 0.2478;
%         disp('HERE!!!! 1 CHECKED AND CORRECTED')
%     end
%     idx4 = find(strcmp({cfg.sensors.sn},'60281'));
%     if ~isempty(idx4)
%          data{idx4}.dn(1831704:end) = data{idx4}.dn(1831704:end) - 0.2372;
%          disp('HERE!!!! 2')
%         plot(data{idx4}.dn)
%     end
% end
end