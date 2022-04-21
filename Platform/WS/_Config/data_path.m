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
        
        ship = matfile(['../../../../data/processed/ShipDas/' Prefix '_ShipDas_Processed.mat']); %%% ship time/location
        
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
        CTD_PROC_Path = [datapath Process_Mode '/RAW_mat_File/']; % saving single profile       
        CTD_PROC_final_Path = [Processed_Path Process_Mode '/']; % combine
        
        ship = matfile(['../../../../data/processed/ShipDas/' Prefix '_ShipDas_Processed.mat']); %%% ship time/location
        
        
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
        
        ADCP_PROC_final_Path = '../../../../data/processed/ADCP/'; % combine
        
        ship = matfile(['../../../../data/processed/ShipDas/' Prefix '_ShipDas_Processed.mat']); %%% ship time/location
       
    case 'PostProcessing'
        %% vmp processing
        Bundle_Path = '../../Processed Bundle/';
        
        ship = load(['../../Processed Bundle/' basename '_shipdata.mat']); %%% ship time/location
        
        [~,temp] = latlon2xy(ship.das.lat,ship.das.lon,lat_c,lon_c);
        
        ship.das.dist = [0;cumsum(temp)];
        
        load([Bundle_Path '/' basename '_combo.mat'])
        
    case 'SectionProcessing'
        
end

