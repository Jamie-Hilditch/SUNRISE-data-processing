clear all, close all

% Locate survey metadata files
dir_in = '../../data/intermediate/survey_metadata';
file_out = '../../section_status.xlsx';
survey_files = dir(fullfile(dir_in,'survey*.csv'));

% Vessels are named as follows in these files:
vessels = {'Pelican', 'Walton_Smith', 'Aries', 'Polly'};

%% Create a spreadsheet file with a single page for each vessel

% 1) Load all survey section data, merge into single table
sections = [];
for i = 1:length(survey_files)
    file_in = fullfile(survey_files(i).folder,survey_files(i).name);
    tmp = readtable(file_in); % Read survey metadata
    tmp.survey = i*ones(height(tmp),1); % Add column for survey number

    % Store or concatenate survey metadata
    if isempty(sections)
        sections = tmp;
    else
        sections = cat(1,sections,tmp);
    end
end

% Rename "n" column to "section"
sections = renamevars(sections,"n","section");

% 2) Write a sheet for each vessel
% Save the following variables
varnames = {'survey','section',...
            'start_utc','start_lat','start_lon',...
            'end_utc','end_lat','end_lon'};

for i = 1:length(vessels)
    % Limit sections to vessel
    vessel_sections = strcmp(sections.vessel,vessels{i});
    vessel_tbl = sections(vessel_sections,varnames);

    % Add empty columns for instrument status
    vessel_tbl.vmp_status = repmat({''},height(vessel_tbl),1);
    vessel_tbl.adcp_status = repmat({''},height(vessel_tbl),1);
    vessel_tbl.tchain_status = repmat({''},height(vessel_tbl),1);
    vessel_tbl.other_status = repmat({''},height(vessel_tbl),1);

    % Write to sheet
    writetable(vessel_tbl,file_out,'Sheet',vessels{i});
end
