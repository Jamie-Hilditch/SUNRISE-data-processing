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

% Also turn off a warning about inefficient partial loading if not using v7.3
warning('off','MATLAB:MatFile:OlderFormat');

% Location of intermediate processed data
proc_dir = fullfile(master_config.data_directory,'processed');

% Location to save section files
dir_out = fullfile(master_config.data_directory,'sections');

% Get survey metadata
surveys = dir(fullfile(proc_dir,'survey_metadata','survey_*_sections.csv'));

% add instrument classes to path
addpath(genpath('Instruments'));


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
PE_300 = ADCP_Ship_Combo('ADCP_PE_wh300',PE_ADCP_file,'wh300');
% setup Pelican hydro
PE_Hydro = Hydro_Combo('HYDRO_Pelican',fullfile(proc_dir,'combined_hydro','SUNRISE2021_PE_hydro_combo.mat'));
% setup PE Tchain
PE_TChain = Tchain('TCHAIN_Pelican',fullfile(Tchain_directory,'Pelican'));

% create the Pelican Vessel class and put it in a vessels structure
vessels.Pelican = Vessel(PE_1200,PE_600,PE_600_no4,PE_300,PE_Hydro,PE_TChain);
fprintf("Setup Pelican Instruments\n")

% setup the Walton Smith Instruments
% setup the Walton Smith ADCPs
WS_ADCP_file = fullfile(proc_dir,'adcp_ship','SUNRISE2021_WS_adcp.mat');
WS_1200 = ADCP_Ship_Combo('ADCP_WS_wh1200',WS_ADCP_file,'wh1200');
WS_600 = ADCP_Ship_Combo('ADCP_WS_wh600',WS_ADCP_file,'wh600');
WS_600_no2 = ADCP_Ship_Combo('ADCP_WS_wh600_no2',WS_ADCP_file,'wh600_no_beam2');
% setup Walton Smith VMP
WS_VMP = VMP_Combo('VMP_Walton_Smith',fullfile(proc_dir,'VMP','SUNRISE2021_WS_combo.mat'));
% setup WS_Tchain
WS_TChain = Tchain('TCHAIN_Walton_Smith',fullfile(Tchain_directory,'Walton_Smith'));

% create WS Vessel Class
vessels.Walton_Smith = Vessel(WS_1200,WS_600,WS_600_no2,WS_VMP,WS_TChain);
fprintf("Setup Walton Smith Instruments\n")

% setup the rhib instuments
% adcps first
ADCP_Polly = ADCP_Rhib('ADCP_Polly',fullfile(proc_dir,'adcp_rhib','polly'),rhib_adcp_variables);
ADCP_Aries = ADCP_Rhib('ADCP_Aries',fullfile(proc_dir,'adcp_rhib','aries'),rhib_adcp_variables);
% rhib Tchains
Polly_TChain = Tchain('TCHAIN_Polly',fullfile(Tchain_directory,'Polly'));
Aries_TChain = Tchain('TCHAIN_Aries',fullfile(Tchain_directory,'Aries'));

% create Vessel Classes for the two rhibs
vessels.Polly = Vessel(ADCP_Polly,Polly_TChain);
fprintf("Setup Polly Instruments\n")
vessels.Aries = Vessel(ADCP_Aries,Aries_TChain);
fprintf("Setup Aries Instruments\n")

% create a string array of vessel names
% NB this code assumes that these names are used in the sections metadata table
% if not then we need a structure mapping from these names to the those in the table
vessel_names = string(fieldnames(vessels))'

% now we loop over surveys
for s = 1:length(surveys)

  fprintf("├──Beginning survey %d\n",s)

  % define the filepath
  file_name = sprintf('SUNRISE_2021_survey_%02d_section_%02d.mat',s,n);
  file_path = fullfile(survey_directory,file_name);

  % delete the old file if it exists
  if isfile(file_path); del(file_path); end;

  % load in survey sections metadata
  sections = readtable(fullfile(surveys(s).folder,surveys(s).name));

  % get section numbers
  section_numbers = unique(sections.n);

  % create a summary table
  summary = table;

  % define the survey directory for section files
  survey_directory = fullfile(dir_out,sprintf('survey_%02d',s));
  if ~exist(survey_directory,'dir'); mkdir(survey_directory); end

  % now loop through different sections
  for n = section_numbers'

    fprintf('  ├── Section %d\n',n);

    % get the relevant rows of the table
    this_section = sections(sections.n == n,:);

    % now loop over the vessels
    for vname in vessel_names
      % get the section start and end times are this vessel
      % these arrays are empty if the vessel is not found for that section
      start_stop = this_section{find(ismember(vname,this_section.vessel)),["start_utc","end_utc"]}

      % if empty we are done with this vessel
      if isempty(start_stop); continue; end;

      % get section data for this vessel
      section_data = vessels.(vname).get_all_data(start_stop(1),start_stop(2));

      % add this data to the file
      save(file_path,'-struct','section_data','-v7.3','-append')

    end

    fprintf('    ├── Saved section to %s\n',file_name);

    % add entry into summary file
    % Record section number
    summary(n,{'Section Number'}) = {n};
    for field = fieldnames(section_data)
      summary(n,field) = {1};
    end


  end

  % save summary file
  summary_file_name = sprintf('SUNRISE_2021_survey_%02d_summary.csv',s);
  writetable(summary,fullfile(survey_directory,summary_file_name));
  fprintf('  ├── Saved survey summary to %s\n',summary_file_name);

end
