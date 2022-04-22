%% The Path file for all files needed for vmp-processing and the

%% name setup
Project = 'SUNRISE2021';
Platform = 'PE';
Prefix = [Project '_' Platform];
datapath = ['../../../../data/Platform/' Platform '/'];
Processed_Path = '../../../../data/Processed/';

%% lat/lon center
lon_c = -92.5;
lat_c =  28.5;
addpath(genpath('../../../toolbox/misc'))

if exist('Process_Mode','var')
    %% crucial files
    switch Process_Mode
        case 'VMP'
            %% vmp processing
            %%% path
            VMP_RAWP_Path = [datapath Process_Mode '/RAW_P_File/']; % raw P file
            VMP_RAWM_Path = [datapath Process_Mode '/RAW_mat_File/']; % raw mat file
            VMP_PROC_Path = [datapath Process_Mode '/Proc_profile_matfile/']; % saving single profile
            VMP_PROC_P_Combine_Path = [datapath Process_Mode '/Proc_combine_matfile/']; % Partially combine
            VMP_PROC_final_Path = [Processed_Path Process_Mode '/']; % combine
            
            ship = matfile([Processed_Path 'ShipDas/' Prefix '_ShipDas_Processed.mat']); %%% ship time/location
            
            %%% loading toolbox
            % add ODAS library to path; use odas V4.4 (please find the newest library
            % in "useful files" in the google drive)
            % ex: addpath('../odas_v4.01/')
            addpath(genpath('../../../toolbox/instrument/odas'))
            addpath(genpath('../../../toolbox/general/gsw'))
            
            
        case 'CTD'
            %% CTD processing
            %%% CTD Path
            CTD_RAWR_Path = [datapath Process_Mode '/RAW_RSK_File/']; % raw rsk file
            CTD_RAWM_Path = [datapath Process_Mode '/RAW_mat_File/']; % raw mat file
            CTD_PROC_Path = [datapath Process_Mode '/Proc_mat_File/']; % saving single profile
            CTD_PROC_final_Path = [Processed_Path Process_Mode '/']; % combine
            
            ship = matfile([Processed_Path 'ShipDas/' Prefix '_ShipDas_Processed.mat']); %%% ship time/location
            
            %%% loading toolbox
            addpath(genpath('../../../toolbox/instrument/rbr-rsktools'))
            addpath(genpath('../../../toolbox/general/gsw'))
            
        case 'ShipDas'
            %% Das processing
            %%% Das Path
            DAS_RAW_Path = [datapath Process_Mode '/']; % raw das file
            DAS_PROC_final_Path = [Processed_Path Process_Mode '/']; % combine
            
        case 'ADCP_UHDAS'
            ADCP_Project_name = 'PE21_24_merge';
            ADCP_name = {'wh1200','wh600','wh600_nobeam4','wh300'};
            ADCP_PROC_Path = [datapath Process_Mode '/' ADCP_Project_name '/proc/'];
            
            ADCP_PROC_final_Path = [Processed_Path Process_Mode '/']; % combine
            
            ship = matfile([Processed_Path 'ShipDas/' Prefix '_ShipDas_Processed.mat']); %%% ship time/location
            
        case 'Tchain'
            TCn_DATA_Path = [datapath Process_Mode '/'];
            TCn_PROC_final_Path = [Processed_Path Process_Mode '/' Platform '/']; % combine
            
            TCn_GPS_Path = [Processed_Path 'ShipDas/' Prefix '_ShipDas_Processed.mat']; %%% ship time/location
            
            %%% loading toolbox
            addpath(genpath('../../../toolbox/instrument/rbr-rsktools'))
            addpath(genpath('../../../toolbox/general/gsw'))

            %%% loading Tchain Config and hook
            addpath(genpath('./tchain_toolbox'))
            
        case 'HydroCombo'
            Hydro_DATA_Path = [Processed_Path Process_Mode '/'];
            
            vmp_combo = load([Processed_Path  'VMP/' Prefix '_VMP_Processed.mat']);
            ctd_combo = load([Processed_Path  'CTD/' Prefix '_CTD_Processed.mat']);
            
        case 'PostProcessing'
            %% vmp processing
            Bundle_Path = '../../Processed Bundle/';
            
            ship = load(['../../Processed Bundle/' basename '_shipdata.mat']); %%% ship time/location
            
            [~,temp] = latlon2xy(ship.das.lat,ship.das.lon,lat_c,lon_c);
            
            ship.das.dist = [0;cumsum(temp)];
            
            load([Bundle_Path '/' basename '_combo.mat'])
    end
end

