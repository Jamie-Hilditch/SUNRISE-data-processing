% file to generate standard CTD data
% can be thought of as different sub-functions
% B zheng
% Dec. 21, 2020
% - % - % - % - % - % - % - % - % - % - % - % - % - % - %- %
%% convert .rsk file into matlab file
clear all
WWmeta.rbrpath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW2/CTD/';
WWmeta.rbrfile=dir(fullfile(WWmeta.rbrpath,'*.rsk'));
WWmeta.name_rbr = 'SR_WW2_D1';
WWmeta.matpath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW2/CTD/';
WWmeta.propath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW/WW2/CTD/';
WWmeta.gridpath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/WW//WW2/CTD/';
WWmeta.lat = 29.59;
WWmeta.lon = -90.7;

WW_rskread(WWmeta);

%% read matlab files
WWmatread(WWmeta);

%% separate into profiles and then grid
WWprofile(WWmeta,15);  % 30 is subject to change

%% grid the CTD product
WWgrid(WWmeta,0.25);

%% sample plot - t/s/sig/chla
para.tscale = [28, 29];
para.sscale = [30, 37];
para.dscale = [20 24];
para.cscale = [0 200];

ctdplot(WWmeta,para)







