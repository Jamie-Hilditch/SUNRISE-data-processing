%% config_SUNRISE.m
% Usage: Called from get_config('SUNRISE') in BowChain_master
% Description: Creates a basic deployment configuration structure for all
%              BowChain/Tchain deployments.
% Inputs: none
% Outputs: config structure
%
% Author: Dylan Winters (dylan.winters@oregonstate.edu)
% Created: 2021-06-20

function config = config_SUNRISE()

%% Set some global default options for deployments
defaults = struct();

% This is the simplest model for initial processing -- just assume that the
% chain is vertical and compressing like a spring. Get the depth coordinate of
% each sensor by interpolating between pressure sensors.

        defaults.chain_model = 'cm_straight';
%   defaults.chain_model = 'cm_segmented';

% We should do a "dunk" to calibrate sensor clocks eventually.
% defaults.time_offset_method = 'cohere';
% defaults.cohere_interval = [dunk_start_time, dunk_end_time];

% Set this to true to force re-parsing of raw data even if .mat files already
% exist. Useful in case we find mistakes in parsing functions and need to
% correct them.
defaults.raw2mat = false;

% Use the earliest and latest timestamp of all sensors if the start/end time
% aren't specified. Just in case we're missing a section_start_end_time.csv.
defaults.dn_range = [-inf inf];

% Grid-averaging settings
defaults.bin_method = 'none';
defaults.bin_dt = 10;
defaults.bin_dz = 1;
defaults.bin_zlim = [-20 0];

% Modify user_directories.m for your data location. For me this returns
% '/home/data/SUNRISE/Tchain', but everyone's local copy will probably be
% different.
dir_raw = fullfile(user_directories('SUNRISE'),'raw');

%% Create deployment-specific configurations
% This is where having a consistent file structure does 90% of the work for us!
% The expected format should look like this:
%
% └── Tchain
%     └── raw
%         └── Aries
%             └── deploy_20210618
%                 ├── 060088_20210618_2140.rsk
%                 ├── 077416_20210618_1644.rsk
%                 ├── 077561_20210618_1649.rsk
%                 ├── 077565_20210618_1647.rsk
%                 ├── 077566_20210618_2148.rsk
%                 ├── 077568_20210618_2141.rsk
%                 ├── 101179_20210618_2145.rsk
%                 ├── 101180_20210618_2136.rsk
%                 ├── instrument_depths.csv
%                 ├── README.txt
%                 └── section_start_end_time.csv

vessel_raw_folders = dir(fullfile(dir_raw));
vessel_names = setdiff({vessel_raw_folders.name},{'.','..'});

ndep = 0;
for v = 1:length(vessel_names)
    deployments = dir(fullfile(dir_raw,vessel_names{v},'deploy*'));
    for i = 1:length(deployments)
        ndep = ndep + 1;
        config(ndep).name = deployments(i).name;
        config(ndep).vessel = vessel_names{v};
        
        % Read the sensors.csv file for instrument deployment positions
        t = readtable(fullfile(deployments(i).folder,deployments(i).name,'instrument_depths.csv'));
        config(ndep).sensor_sn = num2cell(t.serialnum);
        config(ndep).sensor_pos = t.depth_m_;
        
        % Try to read start & end time
        try
            t = readtable(fullfile(deployments(i).folder,deployments(i).name,'section_start_end_time.csv'));
            config(ndep).dn_range = datenum([t.start_time t.end_time]);
        catch err
            % Default to full sensor time range if this fails
            config(ndep).dn_range = [-inf inf];
        end
        
        % Set raw data directory
        config(ndep).dir_raw = fullfile(deployments(i).folder,deployments(i).name);
        
        % Set processed data directory
        config(ndep).dir_proc = strrep(config(ndep).dir_raw,'/raw/','/processed/')
        
        
    end
end
config = fill_defaults(config,defaults);

%% We can specify any deployment- or vessel-specific settings here
% There are a bunch of possible ways to do this, below is one example:
for i = 1:length(config)
    switch config(i).vessel
        
        % Configure Walton Smith deployments
        case 'Walton_Smith'
            switch config(i).name
                case 'deploy_20210621'
                    config(i).zero_pressure_interval = datenum([2021 06 21 21 04 40; 2021 06 21 23 21 20]);
                    config(i).file_gps = '/home/dw/Data/SUNRISE/ws_gps.mat';
            end
            
            % Configure Pelican deployments
        case 'Pelican'
            switch config(i).name
                case 'deploy_20210621_01'
                    config(i).zero_pressure_interval = datenum([2021 06 21 14 04 40; 2021 06 21 14 15 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                case 'deploy_20210623'
                    config(i).zero_pressure_interval = datenum([2021 06 24 17 11 40; 2021 06 24 17 53 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                case 'deploy_20210626'
                     config(i).zero_pressure_interval = datenum([2021 06 26 03 22 00; 2021 06 26 03 57 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                    case 'deploy_20210627'
                     config(i).zero_pressure_interval = datenum([2021 06 27 18 24 40; 2021 06 27 19 23 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                    
                    case 'deploy_20210628'
                     config(i).zero_pressure_interval = datenum([2021 06 29 00 00 49; 2021 06 29 00 38 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                    case 'deploy_20210629'
                     config(i).zero_pressure_interval = datenum([2021 06 29 18 12 00; 2021 06 29 18 19 30]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                    case 'deploy_20210630'
                     config(i).zero_pressure_interval = datenum([2021 07 01 01 00 49; 2021 07 01 1 04 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                    case 'deploy_20210701_A'
                     config(i).zero_pressure_interval = datenum([2021 07 02 00 00 00; 2021 07 02 01 04 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                    case 'deploy_20210701_B'
                      config(i).zero_pressure_interval = datenum([2021 07 02 03 26 49; 2021 07 02 03 30 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                    case 'deploy_20210701_C'
                     config(i).zero_pressure_interval = datenum([2021 07 02 09 42 24; 2021 07 02 09 46 48]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                    case 'deploy_20210701_D'
                     config(i).zero_pressure_interval = datenum([2021 07 02 22 48 29; 2021 07 02 22 52 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    
                    case 'deploy_20210701_E'
                     config(i).zero_pressure_interval = datenum([2021 07 03 00 51 00; 2021 07 03 00 54 59]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
                    

                       case 'deploy_20210706'
                      config(i).zero_pressure_interval = datenum([2021 07 04 02 02 00; 2021 07 04 02 52 00]);
                    % Set lat and lon
                    config(i).file_gps= '/Users/alesanchez-rios/Documents/SUNRISE CRUISE/Ale-files_Oncruise/Ship Data/SUNRISE2021_PE_shipdata.mat';
            end
    end
end


%% End of config_SUNRISE
% After running
%   >> tchain = BowChain_master('SUNRISE','Aries','deploy_20210618')
% you should have a file structure like this:
%
% └── Tchain
%     ├── processed
%     │   └── Aries
%     │       └── deploy_20210618
%     │           ├── 060088_20210618_2140.mat
%     │           ├── 077416_20210618_1644.mat
%     │           ├── 077561_20210618_1649.mat
%     │           ├── 077565_20210618_1647.mat
%     │           ├── 077566_20210618_2148.mat
%     │           ├── 077568_20210618_2141.mat
%     │           ├── 101179_20210618_2145.mat
%     │           └── 101180_20210618_2136.mat
%     └── raw
%         └── Aries
%             └── deploy_20210618
%                 ├── 060088_20210618_2140.rsk
%                 ├── 077416_20210618_1644.rsk
%                 ├── 077561_20210618_1649.rsk
%                 ├── 077565_20210618_1647.rsk
%                 ├── 077566_20210618_2148.rsk
%                 ├── 077568_20210618_2141.rsk
%                 ├── 101179_20210618_2145.rsk
%                 ├── 101180_20210618_2136.rsk
%                 ├── instrument_depths.csv
%                 ├── README.txt
%                 └── section_start_end_time.csv

% And the output, tchain, looks like this:
%
% >> tchain
% tchain =
%   struct with fields:
%       dn: [1x29111 double]
%        t: [8x29111 double]
%        p: [8x29111 double]
%        s: [8x29111 double]
%        x: [8x29111 double]
%        z: [8x29111 double]
%      pos: [8x1 double]
%     info: [1x1 struct]

% You can run BowChain_master('SUNRISE',vessel_name,deploy_name) like this to
% process an individual deployment, or BowChain_master('SUNRISE') to process all
% deployments.
