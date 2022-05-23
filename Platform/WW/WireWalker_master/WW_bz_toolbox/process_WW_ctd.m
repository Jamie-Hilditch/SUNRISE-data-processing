% file to generate standard CTD data
% can be thought of as different sub-functions
% B zheng
% Dec. 21, 2020
% - % - % - % - % - % - % - % - % - % - % - % - % - % - %- %
%% convert .rsk file into matlab file
clear all
WWmeta.rbrpath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/';
WWmeta.rbrfile=dir(fullfile(WWmeta.rbrpath,'*.rsk'));
WWmeta.name_rbr = 'test_rskread';
WWmeta.matpath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/';
WWmeta.propath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/';
WWmeta.gridpath = '/Users/warrbob/Desktop/UCSD/Research/SUNRISE/data/';
WWmeta.lat = 32.86;
WWmeta.lon = -117.27;

WW_rskread(WWmeta);

%% read matlab files
WWmatread(WWmeta);

%% separate into profiles and then grid
WWprofile(WWmeta,30);

%% grid the CTD product
WWgrid(WWmeta,0.25);

%% sample plot - t/s/sig/chla
para.tscale = [9, 12];
para.sscale = [33.7, 34.1];
para.dscale = [25.5 26.5];
para.cscale = [0 2500];

ctdplot(WWmeta,para)







