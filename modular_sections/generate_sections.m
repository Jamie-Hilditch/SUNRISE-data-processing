%% generate_sections.m
% Usage: generate_sections.m
% Description: Generate a structure with data organized by surveys/sections
% defined in the survey_metadata folder.
% Data is loaded from intermediate processed data sources and organized
% into survey/section-level files using consistent logic and section defintions.

clear all, close all

% Add yaml (and other utility) code to path
addpath(genpath('../utils'))

% define filepaths
master_config = yaml.loadFile('../master_config.yml');

% Turn off a warning when filling tables, the default fill is 0 which is exactly what we want
warning('off','MATLAB:table:RowsAddedExistingVars');

% Location of intermediate processed data
proc_dir = fullfile(master_config.data_directory,'processed');

% Location to save section files
dir_out = fullfile(master_config.data_directory,'sections');

% Get survey metadata
surveys = dir(fullfile(proc_dir,'survey_metadata','survey_*_sections.csv'));

% import instrument classes
import Instruments.ADCP_Ship_Combo
import Instruments.ADCP_Rhib
import Instruments.Hydro_Combo
import Instruments.VMP_Combo
import Instruments.Tchain

% set up the instruments individually in structures by vessel
% We can control exactly which instruments and which variables to save to the section files

% define a couple of variables for convenience
Tchain_directory = fullfile(proc_dir,'tchain');
% the rhib files could do with tidying, all the data is saved twice
% we could only save a subset of the these variables because these files are large
rhib_adcp_variables = ["u","v","w","depth","z","dn","config","ahrs_gyro","computed_heading",...
    "vessel_u","vessel_v","time","heading","pitch","roll","nuc_time","bt_range",...
    "bt_range","bt_vel","bt_time","vessel_heading","vessel_lat","vessel_lon",...
    "vel","echo_intens","corr","ahrs_rotm","processing"];

% setup the Pelican instruments
% setup the Pelican ADCPs
PE_ADCP_file = fullfile(proc_dir,'adcp_ship','SUNRISE2021_PE_ADCP.mat');
PE_1200 = ADCP_Ship_Combo('ADCP_PE_wh1200',PE_ADCP_file,'wh1200');
PE_600 = ADCP_Ship_Combo('ADCP_PE_wh600',PE_ADCP_file,'wh600');
PE_600_no4 = ADCP_Ship_Combo('ADCP_PE_wh600_no4',PE_ADCP_file,'wh600_no4');
PE_300 = ADCP_Ship_Combo('ADCP_PE_wh1200',PE_ADCP_file,'wh300');
% setup Pelican hydro
PE_Hydro = Hydro_Combo('HYDRO_Pelican',fullfile(proc_dir,'combined_hydro','SUNRISE2021_PE_hydro_combo.mat'));
% setup PE Tchain
PE_TChain = Tchain('TCHAIN_Pelican',fullfile(Tchain_directory,'Pelican'));

% create the Pelican Vessel class
Pelican = Vessel(PE_1200,PE_600,PE_600_no4,PE_300,PE_Hydro,PE_Tchain);

% setup the Walton Smith Instruments
% setup the Walton Smith ADCPs
WS_ADCP_file = fullfile(proc_dir,'adcp_ship','SUNRISE2021_WS_ADCP.mat');
WS_1200 = ADCP_Ship_Combo('ADCP_WS_wh1200',WS_ADCP_file,'wh1200');
WS_600 = ADCP_Ship_Combo('ADCP_WS_wh600',WS_ADCP_file,'wh600');
WS_600_no2 = ADCP_Ship_Combo('ADCP_WS_wh600_no2',WS_ADCP_file,'wh600_no_beam2');
% setup Walton Smith VMP
WS_VMP = VMP_Combo('VMP_Walton_Smith',fullfile(proc_dir,'VMP','SUNRISE2021_WS_combo.mat'));
% setup WS_Tchain
WS_TChain = Tchain('TCHAIN_Walton_Smith',fullfile(Tchain_directory,'Walton_Smith'));

% create WS Vessel Class
Walton_Smith = Vessel(WS_1200,WS_600,WS_600_no2,WS_VMP,WS_TChain);

% setup the rhib instuments
% adcps first
ADCP_Polly = ADCP_Rhib('ADCP_Polly',fullfile(proc_dir,'adcp_rhib','polly'),rhib_adcp_variables);
ADCP_Aries = ADCP_Rhib('ADCP_Aries',fullfile(proc_dir,'adcp_rhib','aries'),rhib_adcp_variables);
% rhib Tchains
Polly_TChain = Tchain('TCHAIN_Polly',fullfile(Tchain_directory,'Polly'));
Aries_TChain = Tchain('TCHAIN_Aries',fullfile(Tchain_directory,'Aries'));

% create Vessel Classes for the two rhibs
Polly = Vessel(ADCP_Polly,Polly_TChain);
Aries = Vessel(ADCP_Aries,Aries_TChain);

% now we loop over surveys
for s = 1:length(surveys)

  % load in survey sections metadata
  sections = readtable(fullfile(surveys(s).folder,surveys(s).name));

  % get section numbers
  section_numbers = unique(sections.('n'));

  % create a summary table
  summary = table;

  % now loop through different sections
  for n = section_numbers

    % get the relevant rows of the table
    this_section = sections(sections.n == n,:);

    % get the start and end times for each vessel
    % these arrays are empty if the vessel is not found for that section
    PE_start_stop = this_section{find(ismember('Pelican',this_section.vessel)),["start_utc","end_utc"]};
    WS_start_stop = this_section{find(ismember('Walton_Smith',this_section.vessel)),["start_utc","end_utc"]};
    Polly_start_stop = this_section{find(ismember('Polly',this_section.vessel)),["start_utc","end_utc"]};
    Aries_start_stop = this_section{find(ismember('Aries',this_section.vessel)),["start_utc","end_utc"]};

    Pelican_data = Pelican.get_all_data(PE_start_stop;
    Walton_Smith_data = Walton_Smith.get_all_data(WS_start_stop);
    Polly_data = Polly.get_all_data(Polly_start_stop);
    Aries_data = Aries.get_all_data(Aries_start_stop);

    % concatenate structures
    section_data = catstruct(Pelican_data,Walton_Smith_data,Polly_data,Aries_data);

    % save file
    directory = fullfile(dir_out,sprintf('survey_%02d'));
    if ~exist(directory,'dir'); mkdir(directory); end
    file_name = sprintf('SUNRISE_2021_survey_%02d_section_%02d.mat',s,i);
    save(fullfile(directory,filename),'-struct','section_data','-v7.3')

    % add entry into summary file
    % Record section number
    summary(n,{'Section Number'}) = {n};
    for field = fieldnames(section_data)
      summary(n,field) = 1;
    end

  end

  % save summary file
  directory = fullfile(dir_out,sprintf('survey_%02d'));
  summary_file_name = sprintf('SUNRISE_2021_survey_%02d_summary.csv',s);
  writetable(summary,fullfile(directory,summary_file_name));

end



















%
