%
%% The Path file for all files needed for vmp-processing and the
%%% post-processing.

%% name setup

basename='SUNRISE2021_PE';
%basename='SUNRISE2021_WS';

%% lat/lon center

lon_c = -92.5;
lat_c =  28.5;

%% crucial files

switch dictionary
    case 'vmp'
        %% vmp processing
        VMP_RAWP_Path = '../RAW_P_File/'; % raw P file
        VMP_RAWM_Path = '../RAW_mat_File/'; % raw mat file
        VMP_PROC_Path = '../Proc_profile_matfile/'; % saving single profile
        VMP_PROC_final_Path = '../Proc_combine/'; % combine
        Bundle_Path = '../../''Processed Bundle''/';
        
        ship = load([VMP_PROC_final_Path '/' basename '_shipdata.mat']); %%% ship time/location
        
        [~,temp] = latlon2xy(ship.das.lat,ship.das.lon,lat_c,lon_c);
        
        ship.das.dist = [0;cumsum(temp)];
        
        %%% loading toolbox
        % add ODAS library to path; use odas V4.4 (please find the newest library
        % in "useful files" in the google drive)
        % ex: addpath('../odas_v4.01/')

    case 'ctd'
        %% vmp processing
        CTD_RAWR_Path = '../RAW_RSK_File/'; % raw P file
        CTD_RAWM_Path = '../RAW_mat_File/'; % raw mat file
        CTD_PROC_Path = '../Proc_profile_matfile/'; % saving single profile
        CTD_PROC_final_Path = '../Proc_combine/'; % combine
        Bundle_Path = '../../''Processed Bundle''/';
        
        ship = load([CTD_PROC_final_Path '/' basename '_shipdata.mat']); %%% ship time/location
        
        [~,temp] = latlon2xy(ship.das.lat,ship.das.lon,lat_c,lon_c);
        
        ship.das.dist = [0;cumsum(temp)];
        
    case 'PostProcessing'
        %% vmp processing
        Bundle_Path = '../../Processed Bundle/';
        
        ship = load(['../../Processed Bundle/' basename '_shipdata.mat']); %%% ship time/location
        
        [~,temp] = latlon2xy(ship.das.lat,ship.das.lon,lat_c,lon_c);
        
        ship.das.dist = [0;cumsum(temp)];
        
        load([Bundle_Path '/' basename '_combo.mat'])

    case 'SectionProcessing'
        
end

