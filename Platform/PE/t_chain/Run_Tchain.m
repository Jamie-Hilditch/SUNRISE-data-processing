%% initial setup
clear
close all

%%% setting path and pre-loading files
addpath('../_Config')
Process_Mode = 'Tchain';
data_path %% all data path and library

%%% 

Deployment_name = {'all'}; 

% 'all' for all deployments
% [] for the latest deployment
% Or Specific Deployment Names you want to rerun

return
%%%
%TCn_GPS_Path is Ship GPS. If there are several chuck of data,TCn_GPS_Path
%will be empty. You need to manually write pathes into the config function.
config = config_SUNRISE(Deployment_name,TCn_DATA_Path,TCn_GPS_Path);




