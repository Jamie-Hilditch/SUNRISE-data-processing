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

%%%
%TCn_GPS_Path is Ship GPS. If there are several chuck of data,TCn_GPS_Path
%will be empty. You need to manually write pathes into the config function.
config = config_SUNRISE2021_PE(Deployment_name,TCn_DATA_Path,TCn_GPS_Path);

% BowChain_master(Prefix)

chain_struct = BowChain_master(config);

for i = 1:length(chain_struct)
    chain_struct_save =   chain_struct(i) ;
    save([TCn_PROC_final_Path Prefix '_Tchain_' chain_struct_save.info.config.name '_Processed.mat'],'-struct','chain_struct_save','-v7.3')
end