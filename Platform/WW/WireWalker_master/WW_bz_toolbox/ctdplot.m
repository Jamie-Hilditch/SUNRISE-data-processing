function ctdplot(WWmeta,para)
% plot a raw ctd figure
% bz, jne 15, 2021


load([WWmeta.gridpath,'CTDgrid.mat']);


figure
subplot('position',[0.1,0.72,0.8,0.195]);
time = RBRgrid.std_profiles.time;
z = -RBRgrid.std_profiles.z;
h = pcolor(time,z,RBRgrid.std_profiles.T);
set(h, 'EdgeColor', 'none');
colormap(jet);
caxis([para.tscale(1) para.tscale(2)])
hclb = colorbar;
set(hclb,'location','eastoutside'); %colorbar location
set(hclb,'Position',[0.92 0.72 0.015 0.195]);
danum = min(time);
set(gca,'xtick',danum:(max(time)-danum)/4:max(time),'tickdir','in');
% datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(time) max(time)])
title(['\fontsize{15}',WWmeta.name_rbr,' Temperature'])
set(gca,{'xticklabel'},{[]});
% hxlb = xlabel('Time (DD HH:MM)');
% set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)


subplot('position',[0.1,0.505,0.8,0.195]);
h = pcolor(time,z,RBRgrid.std_profiles.S);
set(h, 'EdgeColor', 'none');
colormap(jet);
caxis([para.sscale(1) para.sscale(2)])
hclb = colorbar;
set(hclb,'location','eastoutside'); %colorbar location
set(hclb,'Position',[0.92 0.505 0.015 0.195]);
danum = min(time);
set(gca,'xtick',danum:(max(time)-danum)/4:max(time),'tickdir','in');
% datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(time) max(time)])
title(['\fontsize{15}',WWmeta.name_rbr,' Salinity'])
set(gca,{'xticklabel'},{[]});


subplot('position',[0.1,0.29,0.8,0.195]);
h = pcolor(time,z,RBRgrid.std_profiles.sig0);
set(h, 'EdgeColor', 'none');
colormap(jet);
caxis([para.dscale(1) para.dscale(2)])
hclb = colorbar;
set(hclb,'location','eastoutside'); %colorbar location
set(hclb,'Position',[0.92 0.29 0.015 0.195]);
danum = min(time);
set(gca,'xtick',danum:(max(time)-danum)/4:max(time),'tickdir','in');
% datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(time) max(time)])
title(['\fontsize{15}',WWmeta.name_rbr,' Density'])
set(gca,{'xticklabel'},{[]});


subplot('position',[0.1,0.075,0.8,0.195]);
h = pcolorjw(time,z,RBRgrid.std_profiles.chla);
set(h, 'EdgeColor', 'none');
colormap(jet);
caxis([para.cscale(1) para.cscale(2)])
hclb = colorbar;
set(hclb,'location','eastoutside'); %colorbar location
set(hclb,'Position',[0.92 0.075 0.015 0.195]);
danum = min(time);
set(gca,'xtick',danum:(max(time)-danum)/4:max(time),'tickdir','in');
datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(time) max(time)])
title(['\fontsize{15}',WWmeta.name_rbr,' Chla'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)
end