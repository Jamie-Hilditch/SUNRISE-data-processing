clear
close all

init_period = 2*pi/(2*2*pi/86400*sind(28));

[~, pos1] = tight_subplot(8, 3, [2e-2 2e-3], [3.5e-2 3e-2], [1.5e-2 1.5e-2]);

%figure
%set(gcf, 'PaperUnits', 'inches','PaperPosition', [0 0 18 13.5],'units','inches','position',[0 0 18 13.5]);

for sec_i = 53%51:58

    uneven  = matfile(['../../data/2021 SUNRISE@2021/uneven data/sec_' sprintf('%04i',sec_i) '_raw.mat']);
    gridded = matfile(['../../data/2021 SUNRISE@2021/grid data/sec_' sprintf('%04i',sec_i) '_grid.mat']);
    
    load('../202112 Field_Processing_Dev/post bundle/SUNRISE2021_PE_ADCP.mat')
    
    L_ship_timestep = uneven.L_ship_timestep;
    
    wh600_t = uneven.wh600_t;
    up_600 = uneven.up_600;
    uf_600 = uneven.uf_600;
    u_600 = uneven.u_600;
    v_600 = uneven.v_600;   
    
    wh1200_t = uneven.wh1200_t;
    up_1200 = uneven.up_1200;
    uf_1200 = uneven.uf_1200;
    u_1200 = uneven.u_1200;
    v_1200 = uneven.v_1200;     
    
    
    lat = uneven.lat;
    lon = uneven.lon;
    heading = uneven.heading;
    
    x = 2*pi*6378e3*cosd(29)/360;
    y = 2*pi*6378e3/360;
    
    u_geo = diffxy(L_ship_timestep,lon)*x/86400;
    v_geo = diffxy(L_ship_timestep,lat)*y/86400;
    
    %%
    close all
    %%
    figure
    subplot(2,1,1)
    plot(wh600_t,wrapTo360(angle(uf_600(1,:)+1i*up_600(1,:))/(2*pi)*360))
    hold on
    plot(wh600_t,wrapTo360(angle(u_600(1,:)+1i*v_600(1,:))/(2*pi)*360))
    hold off
    title('wh600 (bin=1)')
    ylabel('degree')
    legend('ang(uf+iup)','ang(u+iv)') 
    xlim([min(L_ship_timestep) max(L_ship_timestep)])           
    datetick('x','mm/dd HHMM','keeplimits')
    
    subplot(2,1,2)
    plot(wh600_t,wrapTo360(angle(u_600(1,:)+1i*v_600(1,:))/(2*pi)*360)-wrapTo360(angle(uf_600(1,:)+1i*up_600(1,:))/(2*pi)*360))
    hold on
    plot(wh600_t,wrapTo360(angle(u_600(5,:)+1i*v_600(5,:))/(2*pi)*360)-wrapTo360(angle(uf_600(5,:)+1i*up_600(5,:))/(2*pi)*360))
    plot(wh600_t,wrapTo360(angle(u_600(10,:)+1i*v_600(10,:))/(2*pi)*360)-wrapTo360(angle(uf_600(10,:)+1i*up_600(10,:))/(2*pi)*360))
    plot(L_ship_timestep,movmean(angle(u_geo+1i*v_geo)/(2*pi)*360,60))
    plot(L_ship_timestep,wrapTo180(winddeg_conv(heading,1)-180))
    hold off
    title('heading')
    ylabel('degree') 
    legend('bin=1','bin=5','bin=10','real','GPS Heading')
    xlim([min(L_ship_timestep) max(L_ship_timestep)])           
    datetick('x','mm/dd HHMM','keeplimits')
 
    %%
    figure
    subplot(2,1,1)
    plot(wh1200_t,wrapTo360(angle(uf_1200(1,:)+1i*up_1200(1,:))/(2*pi)*360))
    hold on
    plot(wh1200_t,wrapTo360(angle(u_1200(1,:)+1i*v_1200(1,:))/(2*pi)*360))
    hold off
    ylabel('degree')
    legend('ang(uf+iup)','ang(u+iv)')
    title('wh1200 (bin=1)')
    xlim([min(L_ship_timestep) max(L_ship_timestep)])       
    datetick('x','mm/dd HHMM','keeplimits')
    
    subplot(2,1,2)
    plot(wh1200_t,wrapTo360(angle(u_1200(1,:)+1i*v_1200(1,:))/(2*pi)*360)-wrapTo360(angle(uf_1200(1,:)+1i*up_1200(1,:))/(2*pi)*360))
    hold on
    plot(wh1200_t,wrapTo360(angle(u_1200(5,:)+1i*v_1200(5,:))/(2*pi)*360)-wrapTo360(angle(uf_1200(5,:)+1i*up_1200(5,:))/(2*pi)*360))
    plot(wh1200_t,wrapTo360(angle(u_1200(10,:)+1i*v_1200(10,:))/(2*pi)*360)-wrapTo360(angle(uf_1200(10,:)+1i*up_1200(10,:))/(2*pi)*360))
    plot(L_ship_timestep,movmean(angle(u_geo+1i*v_geo)/(2*pi)*360,60))
    plot(L_ship_timestep,wrapTo180(winddeg_conv(heading,1)-180))    
    hold off   
    title('heading')
    ylabel('degree')    
    legend('bin=1','bin=5','bin=10','real','GPS Heading')
    xlim([min(L_ship_timestep) max(L_ship_timestep)])       
    datetick('x','mm/dd HHMM','keeplimits')

    %%    
    figure
    subplot(3,1,1)    
    plot(wh600.dn,-wh600.umeas_bt)
    hold on
    plot(L_ship_timestep+30/86400,movmean(u_geo,60))
    hold off
    title('u ship velocity')
    legend('bottom track','real')
    xlim([min(L_ship_timestep) max(L_ship_timestep)])
    datetick('x','mm/dd HHMM','keeplimits')
    
    subplot(3,1,2)    
    plot(wh600.dn,-wh600.vmeas_bt)
    hold on
    plot(L_ship_timestep+30/86400,movmean(v_geo,60))
    hold off
    title('v ship velocity')
    legend('bottom track','real')    
    xlim([min(L_ship_timestep) max(L_ship_timestep)])    
    datetick('x','mm/dd HHMM','keeplimits')
    
    subplot(3,1,3)   
    u_temp = interp1(L_ship_timestep+30/86400,movmean(u_geo,60),wh600.dn);
    v_temp = interp1(L_ship_timestep+30/86400,movmean(v_geo,60),wh600.dn);
    
    plot(wh600.dn,wrapTo360(angle(u_temp+1i*v_temp)/(2*pi)*360)-wrapTo360(angle(-wh600.umeas_bt+1i*-wh600.vmeas_bt)/(2*pi)*360))
    xlim([min(L_ship_timestep) max(L_ship_timestep)])   
    title('ang(bt) - ang(real)')
    ylabel('deg')
    datetick('x','mm/dd HHMM','keeplimits')
    
end