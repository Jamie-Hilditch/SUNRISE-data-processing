function plot_result_adcp(WWmeta,Vel)
% code to take a quick look at the result
% Bofu Zheng

cb = cbrewer('div','RdBu',21);
figure
subplot(211)
time = Vel{3};
dz = -Vel{4};
h = pcolor(time,dz,Vel{1});
set(h, 'EdgeColor', 'none');
colormap(flipud(cb));
caxis([-0.2 0.2])
colorbar
danum = min(min(time));
set(gca,'xtick',danum:1:max(max(time)),'tickdir','in');
datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(min(time)) max(max(time))])
title(['\fontsize{15}',WWmeta.name_aqd,' E-W velocity'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)


subplot(212)
h = pcolor(time,dz,Vel{2});
set(h, 'EdgeColor', 'none');
colormap(flipud(cb));
caxis([-0.2 0.2])
colorbar
danum = min(min(time));
set(gca,'xtick',danum:1:max(max(time)),'tickdir','in');
datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(min(time)) max(max(time))])
title(['\fontsize{15}',WWmeta.name_aqd,' N-S velocity'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)
end
