% main code to process WW ADCP and CTD data
% ADCP is in downward looking configuration
%
% Bofu Zheng/ Arnaud Le Boyer/ Drew Lucas
% Nov.12 2020
% boz080@ucsd.edu

%% set path
% here the code is fed with .mat files of the measured velocity. To obtain
% .mat files, user need to run Signature Deployment software first to
% convert .ad2cp files into .mat files

clc; clear all;
WWmeta.aqdpath='/Users/warrbob/Desktop/ECOHAB/converted/'; % path for raw .mat Nortek Signature data
WWmeta.root_script='/Users/warrbob/Desktop/UCSD/Research/WireWalker_master/'; % root for WW_master toolbox
WWmeta.name_aqd=['ECOHAB1']; % Name of the processed data, can be changed according to different cruises
WWmeta.matpath='/Users/warrbob/Desktop/ECOHAB/mat/'; % path to save the processed data - mat
WWmeta.propath='/Users/warrbob/Desktop/ECOHAB/pro/'; % path to the save processed data - profile
WWmeta.propath_rearrange='/Users/warrbob/Desktop/ECOHAB/pro_re/'; % path to the copied profiles - for combining cut-off profiles
WWmeta.gridpath='/Users/warrbob/Desktop/ECOHAB/grid/'; % path to the save estimated velocity field
WWmeta.figpath='/Users/warrbob/Desktop/ECOHAB/fig/'; % (optional) path to save the profiling figure, can be commented out
WWmeta % display what has been entered
cd(WWmeta.root_script) % change directory to the location...
dd0 = dir([WWmeta.aqdpath '*.mat']);

%% set variables
% adjustable variables include:
% 
% NUM_combining_files: number of .mat files to be combined as a group from 
%                      raw output of the ADCP, typically is set to be 20.
%                      if this number is too large, combined file is too 
%                      big to save
% blockdis           : blocking distance, can be found in the Config
%                      structure of the raw mat file
% cellsize           : cell size, can be found in the Config structure of
%                      the raw mat file
% saprate            : sampling rate, in Hz, can be found in the Config
%                      structure of the raw mat file
% boxsize            : vertical range for averaging 
%                      - determining the vertical resolution of the final product
%                      typically is set to be the same as cell size
% z_max              : max depth of the WW profile, positive value
% k                  : determine whether to process downcast data. if k==1,
%                      downcast data will be saved. if k~=1, only upcast
%                      data will be processed and saved.
% 
variables.NUM_combining_files = 20;  % need to specify the number for combining files, the default is: 20
variables.blockdis = 0.1;   % need to specify the blocking distance, the default is: 0.1m
variables.cellsize = .5;     % need to specify the cell size, the default is: 0.25m
variables.saprate = 16;     % need to specify the sampling rate, the default is: 16Hz
variables.boxsize = 0.5;      % is set to be 1m vertically, default is 0.25m
variables.z_max   = 55;    % 500m profile, the default is: 100m
k = 1;                     % 1 or not 1


%% sort files
% this is to make sure raw .mat files are combined in the right order
WWmeta = sort_file(WWmeta)

%% combine separate raw .mat files together and then chunk into profiles
for q = 1:variables.NUM_combining_files:length(dd0)  
    if q+variables.NUM_combining_files-1>length(dd0)  % if length of the index is longer than the number of the raw data, then the ending index is set to be the lenght of the raw data file
        num = length(dd0)-q+1;
    else
        num = variables.NUM_combining_files;
    end
    
    merge_signature(WWmeta,q,num);    % merge separate .mat files from ADCP output and then form a group
    create_profiles(WWmeta,q,num,200,k);  % will chunk ww profiles into upcast and downcast, and save them separately
    disp(['current file location: ',num2str(q),'_',num2str(q+num-1)])  % show where we are
end
disp('identify profiles: finished')


%% combine cut-off profiles
% there may be some profiles (specifically the first profile or last profile in a group)
% with first half in the previous group and second half in the current group. 
% Therefore, we are going to combine the cut-off profiles. 
copyfile([WWmeta.propath,'*.mat'],WWmeta.propath_rearrange)  % copy file from the old folder to the new folder
disp('copying file: finished')
combine_cutoff(WWmeta,variables.NUM_combining_files,k)  % combine_cutoff is performed in the new folder
disp('combining: finished')

%% WWvel analysis
% here the motion correction and box averaging are perfomed
Vel = WWvel_downward(WWmeta,variables);  % function to generate estimated velocity field
Vel{1}(find(abs(Vel{1})>0.3)) = nan;   % optional, remove outliers
Vel{2}(find(abs(Vel{2})>0.3)) = nan;   % optional, remove outliers
ADCP.time = nanmean(Vel{3},1);
ADCP.dz   = Vel{4};
ADCP.velE = Vel{1};
ADCP.velN = Vel{2};

save([WWmeta.gridpath,WWmeta.name_aqd,'_vel.mat'],'ADCP');  % save result

%% take a quick look at the result
plot_result_adcp(WWmeta,Vel)



%% process ctd data

%% change to path of WW_master
WWctd.root_script='/Users/warrbob/Desktop/UCSD/Research/WireWalker_Navo-master/WireWalker_Navo-master/'; % root for WW_master toolbox
WWctd.pathcervello='/Volumes/warrbob3/misobob/DATA/miso3/DATA/CTD/RAW/Cervello/RBR-110050.csv'; % path for raw cervello data
WWctd.pathconcerto='/Volumes/warrbob3/misobob/DATA/miso3/DATA/CTD/RAW/Concerto/RBR-66044.csv';  % path for raw concerto data
WWctd.matpath = '/Volumes/warrbob3/misobob/DATA/miso3/DATA/CTD/mat/';   % path for saving raw .mat data
WWctd.propath = '/Volumes/warrbob3/misobob/DATA/miso3/DATA/CTD/pro/';   % path for saving .mat data with profiles
WWctd.gridpath = '/Volumes/warrbob3/misobob/DATA/miso3/DATA/CTD/grid/'; % path for saving .mat data with grids
WWctd; % display what has been entered
cd(WWmeta.root_script)


%% process telemetry CTD data
% read data from two different folders
disp(['read ' WWctd.pathcervello])
Cervello = mod_telemetry_Cervello_struct(WWctd.pathcervello);
disp(['read ' WWctd.pathconcerto])
Concerto = mod_telemetry_Concerto_struct(WWctd.pathconcerto);
% combine Cervello and Concerto data together (interpolate Cervelle on to
% Concerto)
CTD = mod_merge_Concerto_Cervello(Cervello,Concerto);

% get profiles
% min_depth is set to be a weird number, otherwise, getcastctd
% function will generate uo with no data inside, then function breaks
[up,down,dataup,datadown] = mod_telemetry_getcastctd(CTD,3.1,10);
% get grids
CTDgrid = mod_telemetry_gridctd(dataup,0.25,100); % grid the Wirewalker data to a depth and time grid.
... Depth bins are 0.25 m by default. Time is the center time point of each profile.
    ...Note that this means that time is NOT even spaced in this approach to gridding.


%% save mat file
% file name sp
MatFILENAME = [Concerto.sn,'.mat'];
ProfileFILENAME = ['Profile_',Concerto.sn,'_','.mat'];
GridFILENAME = ['Grid_',Concerto.sn,'_','.mat'];

save([WWctd.matpath,MatFILENAME],'CTD');
save([WWctd.propath,ProfileFILENAME],'dataup');
save([WWctd.gridpath,GridFILENAME],'CTDgrid');

