% put DO data onto WW data set
clear all
WWmeta.rbrpath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/CTD/';
WWmeta.rbrfile=dir(fullfile(WWmeta.rbrpath,'*.rsk'));
WWmeta.name_rbr = 'SR_WW1_D2';
WWmeta.matpath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/CTD/';
WWmeta.propath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/CTD/';
WWmeta.gridpath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/CTD/';
WWmeta.lat = 29.59;
WWmeta.lon = -90.7;

% DO data
data = readtable(['/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/DO/SUNRISE_WW1_D2_ARO-USB_0186_191948_A.csv']);
%time
time = datenum(data{:,1});
DO = data{:,4};   % in mg/L
Temp = data{:,2};

%% load WWdata
load('/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/CTD/CTDall.mat')

%% find unique time for DO sensor
[~,I] = unique(time);
ind = 1:length(time);
ind(I) = [];
time(ind) = time(ind-1)+0.5/86400;

%% interpolation
ind = find(diff(time)>0.015);
CTDall.DO = interp1(time,DO,CTDall.time);
for i = 1:length(ind)
    indnan =  find(CTDall.time>time(ind(i)) & CTDall.time<time(ind(i)+1));
    CTDall.DO(indnan) = nan;
end

save([WWmeta.matpath,'CTDall.mat'],'CTDall');

%%
WWprofile(WWmeta,15);  

%%
WWgrid(WWmeta,0.25);


%% for two hobos
clear all
data = readtable(['/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/uppper hobo/PISCO_2_deploy2.csv']);
% time
time = datenum(data{:,2});
DO = data{:,3};   % in mg/L
Temp = data{:,4};
time = time + datenum(2021,6,30,12,0,0) - time(1);  % adjust time

hobo_up.time = time+5/24;
hobo_up.DO = DO;
hobo_up.temp = Temp;
save('/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/uppper hobo/hobo_up.mat','hobo_up');



%% lower hobo
clear all
data = readtable(['/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/lower hobo/Shearman_2_deploy2.csv']);
% time
time = datenum(data{:,2});
DO = data{:,3};   % in mg/L
Temp = data{:,4};
time = time + datenum(2021,6,30,12,0,0) - time(1);  % adjust time

hobo_low.time = time+5/24;
hobo_low.DO = DO;
hobo_low.temp = Temp;

save('/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/lower hobo/hobo_low.mat','hobo_low');


%% plot to take a look
clear all

load('/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/uppper hobo/hobo_up.mat');
load('/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/lower hobo/hobo_low.mat');
load('/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW1/D2/CTD/CTDgrid.mat');

cbo   = cmocean('oxy');
danum = datenum(2021,6,29);

figure
subplot('Position',[0.1 0.42 0.78 0.5]);
pcolorjw(RBRgrid.std_profiles.time,-RBRgrid.std_profiles.z,RBRgrid.std_profiles.DO)
colormap(cbo)
caxis([6.5 8])
hclb = colorbar;
set(hclb,'location','eastoutside'); %colorbar location
set(hclb,'Position',[0.9 0.42 0.01 0.5],'ytick',3:0.5:8.5,...
    'tickdir','out');
set(gca,'xtick',danum:1:RBRgrid.std_profiles.time(end),'tickdir','in');
set(gca,{'xticklabel'},{[]});
xlim([RBRgrid.std_profiles.time(1) RBRgrid.std_profiles.time(end)])
ylabel('depth (m)')
title('WW1 deploy2 oxygen (mg L^{-1})')


subplot('Position',[0.1 0.1 0.78 0.3]);
plot(RBRgrid.std_profiles.time,nanmean(RBRgrid.std_profiles.DO(end-4:end,:),1),'b.-')
hold on
plot(hobo_up.time,hobo_up.DO,'k.-')
plot(hobo_low.time,hobo_low.DO,'r.-')

legend('WW bottom 1 m averaged','upper hobo','lower hobo','location','southwest')
set(gca,'xtick',danum:1:RBRgrid.std_profiles.time(end),'tickdir','in');
datetick('x','mm/DD','keepticks');
xlim([RBRgrid.std_profiles.time(1) RBRgrid.std_profiles.time(end)])
ylim([6 8.5])
ylabel('oxygen (mg L^{-1})')
xlabel('time (month/day in 2021)')



