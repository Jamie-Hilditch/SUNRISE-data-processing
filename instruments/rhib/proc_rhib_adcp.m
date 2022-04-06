% This is the script I used to parse/process Nortek ADCP, GPS, and IMU data
% from the RHIBs in 2021.
%
% This is a high-level script around some high-level functions that call
% lower-level functions in the ROSE_code directory. I've added some comments
% explaining what code gets called where.
%
% Note: This is a copy of what was on my laptop after the cruise, plus some
% additional comments. I've since made some improvements to a handful of the
% low-level functions, which can be found here:
% - https://github.com/dswinters/ocean-tools
%
% I recommend using this code as a starting point for building a new script
% around the lower-level tools.


clear all, close all
addpath('ROSE_code');

output_dir = ''; % where to save processed .mat files
raw_dir = ''; % where to find deployment folders, i.e. the
              % directory containing these folders:
% depname = 'UBOX01_deploy_20210702_141010';
% depname = 'UBOX01_deploy_20210704_214603';
% depname = 'UBOX01_deploy_20210706_215529';
% depname = 'UBOX02_deploy_20210629_102001';
% depname = 'UBOX02_deploy_20210630_202914';
% depname = 'UBOX02_deploy_20210630_202914';
% depname = 'UBOX02_deploy_20210701_234811';
% depname = 'UBOX02_deploy_20210704_005204';

% These are split into RHIB-specific directories on Drive under:
% NIW_GoM_Processing/2021/data/raw/rhib/Aries
% NIW_GoM_Processing/2021/data/raw/rhib/Polly

% Pick a deployment. I was processing deployments as RHIBs were recovered, so
% doing it like this made sense at the time. For re-processing all the RHIB
% data, definitely turn this whole script into a loop over all deployments
% instead.
depname = 'UBOX02_deploy_20210706_002617';

%% Parse or load raw data
disp(depname)
disp('Getting raw data...')
raw_file = fullfile(proj_dir,'mat',[depname '.mat']);
overwrite_raw = false; % set to true to always parse from raw data
% Otherwise, parse raw data only if no .mat version exists; else load the .mat version.
if ~exist(raw_file,'file') | overwrite_raw;

    % This function parses all of the raw data by calling other functions to
    % parse specific instrument data. See:
    % - rose_parse_deployment.m
    %   - parse_nortek_adcp.m
    %     - nortek_parse_bt.m
    %     - nortek_parse_burst.m
    %   - parse_gps.m
    %   - parse_imu.m
    data_raw = rose_parse_deployment(fullfile(raw_dir,depname));

    save('-v7.3',fullfile(proj_dir,'mat',[depname '.mat']),'-struct','data_raw');
else
    data_raw = load(raw_file);
end
disp('Done!')

%% Process ADCP data
disp('Processing...')
% See ROSE_code/rose_proc_adcp for ADCP processing steps. Improvements to
% processing (better filtering, motion handling, etc.) should be made there.
%
% It relies on inputs specifying how to compute ADCP orientation and vessel
% velocity so that different methods can be substituted when instrumentation
% changes. See:
% - ROSE_code/compute_adcp_orientation.m
% - ROSE_code/compute_vessel_vel.m
% for how these are actually computed based on the methods specified here:
adcp_orientation_method = {'Nortek AHRS + Hemisphere + offset', 45};
vessel_vel_method = 'GPRMC groundspeed and course';
adcp = rose_proc_adcp(data_raw,adcp_orientation_method,vessel_vel_method);

disp('Done!')

%% Save processed data
disp('Saving...')
save('-v7.3',fullfile(proj_dir,'proc',['adcp_' depname '.mat']),'-struct','adcp');
disp('Done!')
