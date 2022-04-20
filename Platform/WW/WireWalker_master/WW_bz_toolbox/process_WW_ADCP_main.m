% main code to process WW ADCP data
% ADCP is in downward looking configuration
%
% Bofu Zheng/ Drew Lucas/ Arnaud Le Boyer/ 
% Oct.11 2021
% boz080@ucsd.edu

%% set path
% here the code is fed with .mat files of the measured velocity. To obtain
% .mat files, user needs to run Signature Deployment software first to
% convert .ad2cp files into .mat files

clc; clear all;
WWmeta.aqdpath='/Users/warrbob/Desktop/ECOHAB/convertd/'; % path for raw .mat Nortek Signature data
WWmeta.root_script='/Users/warrbob/Desktop/UCSD/Research/WireWalker_master/'; % root for WW_ADCP toolbox
WWmeta.name_aqd=['ECOHAB1']; % Name of the processed data, can be changed according to different cruises
WWmeta.matpath='/Users/warrbob/Desktop/ECOHAB/mat/'; % path to save the processed data - mat (combined raw .mat files)
WWmeta.propath='/Users/warrbob/Desktop/ECOHAB/pro/'; % path to the save processed data - profile (separated into up/down casts)
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
%                      big to be saved
% blockdis           : blocking distance, can be found in the Config
%                      structure of the raw .mat file
% cellsize           : cell size, can be found in the Config structure of
%                      the raw .mat file
% saprate            : sampling rate, in Hz, can be found in the Config
%                      structure of the raw .mat file
% boxsize            : vertical range for averaging 
%                      - determining the vertical resolution of the final product
%                      typically is set to be the same as cell size
% z_max              : max depth of the WW profile, positive value
% k                  : determine whether to process downcast data. if k==1,
%                      downcast data will be saved. if k~=1, only upcast
%                      data will be processed and saved.
% thhold             : threshold value to determine whether it is too short
%                      to be a profile
variables.NUM_combining_files = 20;  % need to specify the number for combining files, the default is: 20
variables.blockdis = 0.1;            % need to specify the blocking distance, the default is: 0.1m
variables.cellsize = .5;             % need to specify the cell size, the default is: 0.25m
variables.saprate = 16;              % need to specify the sampling rate, the default is: 16Hz
variables.boxsize = 0.5;             % is set to be 1m vertically, default is 0.25m
variables.z_max   = 55;              % 500m profile, the default is: 100m
variables.k = 0;                     % 1 or not 1
variables.thhold = 200;              % need to specify the number


%% sort files
% this is to make sure raw .mat files are in the right order
WWmeta = sort_file(WWmeta)

%% combine separate raw .mat files together and then chunk into profiles
for q = 1:variables.NUM_combining_files:length(dd0)  
    if q+variables.NUM_combining_files-1>length(dd0)
        num = length(dd0)-q+1;
    else
        num = variables.NUM_combining_files;
    end
    
    merge_signature(WWmeta,q,num);    % merge separate .mat files from ADCP output and then form a group
    create_profiles(WWmeta,q,num,variables.thhold,variables.k);  % will chunk ww profiles into upcast and downcast, and save them separately
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
% here the motion correction and box averaging are performed
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




