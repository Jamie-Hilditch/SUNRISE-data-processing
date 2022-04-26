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
% warning('off','MATLAB:MatFile:OlderFormat');

% Location of intermediate processed data
proc_dir = fullfile(master_config.data_directory,'Processed');

% Location to save section files
dir_out = fullfile(master_config.data_directory,'Sections');

% Get survey metadata
surveys = dir(fullfile(proc_dir,'survey_metadata','survey_*_sections.csv'));

% add instrument classes to path
addpath(genpath('Instruments'));

% ###### SETUP VESSELS AND INSTRUMENTS ###### %

fprintf("Setting up vessels and instruments\n")

% set up the instruments individually in by vessel
% we can control exactly which instruments and which variables to save to the section files

% setup the Pelican instruments
% define Pelican data sources
PE_wh300_file = fullfile(proc_dir,'ADCP','SUNRISE2021_PE_wh300_Processed.mat');
PE_wh600_file = fullfile(proc_dir,'ADCP','SUNRISE2021_PE_wh600_Processed.mat');
PE_wh600_nobeam4_file = fullfile(proc_dir,'ADCP','SUNRISE2021_PE_wh600_nobeam4_Processed.mat');
PE_wh1200_file = fullfile(proc_dir,'ADCP','SUNRISE2021_PE_wh1200_Processed.mat');
PE_Hydro_file = fullfile(proc_dir,'HydroCombo','SUNRISE2021_PE_HydroCombo_Processed.mat');
PE_TChain_directory = fullfile(proc_dir,'Tchain','PE');
PE_FTMET_file = fullfile(proc_dir,'ShipDas','SUNRISE2021_PE_ShipDas_Processed.mat');

% create the Pelican Vessel class and instruments
% store the vessels in a structure
vessels.Pelican = Vessel(ADCP_Ship_Combo('ADCP_PE_wh1200',PE_wh1200_file), ...
                         ADCP_Ship_Combo('ADCP_PE_wh600',PE_wh600_file), ...
                         ADCP_Ship_Combo('ADCP_PE_wh600_no4',PE_wh600_nobeam4_file), ...
                         ADCP_Ship_Combo('ADCP_PE_wh300',PE_wh300_file), ...
                         Hydro_Combo('HYDRO_PE',PE_Hydro_file), ...
                         Tchain('TCHAIN_PE',PE_TChain_directory), ...
                         FTMET_mat('FTMET_PE',PE_FTMET_file));

fprintf("Setup Pelican Instruments\n")

% setup the Walton Smith Instruments
% define the Walton Smith data sources
WS_wh300_file = fullfile(proc_dir,'ADCP','SUNRISE2021_WS_wh300_Processed.mat');
WS_wh600_file = fullfile(proc_dir,'ADCP','SUNRISE2021_WS_wh600_Processed.mat');
WS_wh600_nobeam2_file = fullfile(proc_dir,'ADCP','SUNRISE2021_WS_wh600_nobeam2_Processed.mat');
WS_wh1200_file = fullfile(proc_dir,'ADCP','SUNRISE2021_WS_wh1200_Processed.mat');
WS_Hydro_file = fullfile(proc_dir,'HydroCombo','SUNRISE2021_WS_HydroCombo_Processed.mat');
%WS_TChain_directory = fullfile(proc_dir,'tchain','Walton_Smith');
WS_FTMET_file = fullfile(proc_dir,'ShipDas','SUNRISE2021_WS_ShipDas_Processed.mat');

% create WS Vessel Class and intruments
vessels.Walton_Smith = Vessel(ADCP_Ship_Combo('ADCP_WS_wh1200',WS_wh1200_file), ...
                              ADCP_Ship_Combo('ADCP_WS_wh600',WS_wh600_file), ...
                              ADCP_Ship_Combo('ADCP_WS_wh600_no2',WS_wh600_nobeam2_file), ...
                              Hydro_Combo('HYDRO_WS',WS_Hydro_file), ...
                              FTMET_mat('FTMET_WS',WS_FTMET_file));

fprintf("Setup Walton Smith Instruments\n")

% setup Polly Instruments
%Polly_TChain_directory = fullfile(proc_dir,'tchain','Polly');

% create Polly Vessel Class and intruments
%vessels.Polly = Vessel(Tchain('TCHAIN_Polly',Polly_TChain_directory));

%fprintf("Setup Polly Instruments\n")

% setup Aries Instruments
%Aries_TChain_directory = fullfile(proc_dir,'tchain','Aries');

% create Aries Vessel Class and intruments
%vessels.Aries = Vessel(Tchain('TCHAIN_Aries',Aries_TChain_directory));

%fprintf("Setup Aries Instruments\n")

% create a string array of vessel names
% note that this code assumes that these names are used in the sections metadata table
% if not then we need a structure mapping from these names to the those in the table
vessel_names = string(fieldnames(vessels))';

fprintf("Setup Complete\n\n")

% ###### SETUP COMPLETE ###### %

% ###### NOW GENERATE SOME SECTIONS ###### %

% we loop over surveys
for s = 1:length(surveys)

  fprintf("├──Beginning survey %d\n",s)

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

    % create a cell array to store section data
    section_data = cell.empty;

    % now loop over the vessels
    for vname = vessel_names

      % get the section start and end times are this vessel
      % these arrays are empty if the vessel is not found for that section
      start_stop = this_section{find(ismember(vname,this_section.vessel)),["start_utc","end_utc"]};

      % if empty we are done with this vessel
      if isempty(start_stop); continue; end;

      % get section data for this vessel
      section_data = [section_data,vessels.(vname).get_all_data(start_stop(1),start_stop(2))];

    end

    % concatenate all data into one structure
    section_structure = catstruct(section_data{:});

    % save file
    file_name = sprintf('SUNRISE_2021_survey_%02d_section_%02d.mat',s,n);
    save(fullfile(survey_directory,file_name),'-struct','section_structure','-v7.3')

    fprintf('    ├── Saved section to %s\n',file_name);

    % add a row into the summary table
    summary(n,{'Section Number'}) = {n};
    for field = fieldnames(section_structure)
      summary(n,field) = {1};
    end

  end

  % sort summary table instruments alphabetically
  table_instruments = summary.Properties.VariableNames;
  table_instruments(strcmp(table_instruments,'Section Number')) = [];
  summary = summary(:,['Section Number',sort(table_instruments)]);

  % save the summary table
  summary_file_name = sprintf('SUNRISE_2021_survey_%02d_summary.csv',s);
  writetable(summary,fullfile(survey_directory,summary_file_name));
  fprintf('  ├── Saved survey summary to %s\n',summary_file_name);

end


% rhib_adcp_variables = ["u","v","w","depth","z","dn","config","ahrs_gyro","computed_heading",...
%     "vessel_u","vessel_v","time","heading","pitch","roll","nuc_time","bt_range",...
%     "bt_range","bt_vel","bt_time","vessel_heading","vessel_lat","vessel_lon",...
%     "vel","echo_intens","corr","ahrs_rotm","processing"];
