%% generate_sections.m
% Usage: generate_sections.m
% Description: Generate a structure with data organized by surveys/sections
% defined in the survey_metadata folder.
% Data is loaded from intermediate processed data sources and organized
% into survey/section-level files using consistent logic and section defintions.

clear all, close all

% define filepaths
master_config = yaml.loadFile("master_config.yml")
addpath(genpath(master_config.code_directory))

% Define some file/folder locations
addpath(fullfile('..','general'));

% Location of intermediate processed data
proc_dir = fullfile(master_config.data_directory,'processed');

% Location to save section files
dir_out = fullfile(master_config.data_directory,,'sections');

% Get survey metadata
surveys = dir(fullfile(proc_dir,'survey_metadata','survey_*_sections.csv'));

% Overwrite existing section files?
overwrite_sections = true;

% Define vessel names & short names. This struct is used to assign vessel
% numbers and index by vessel.
vessels = struct(...
    'name',      {'Pelican', 'Walton_Smith', 'Aries', 'Polly'},...
    'shortname', {     'PE',           'WS', 'RHIB1', 'RHIB2'});

% Begin looping over surveys
for s = 1:length(surveys)

    % Initialize output structure
    fprintf('Survey %d\n',s);

    % Get section information from survey
    secs = readtable(fullfile(surveys(s).folder,surveys(s).name));
    sname = surveys(s).name(1:regexp(surveys(s).name,'_sections.csv','start')-1);

    %% Load single-file data sources
    % UHDAS ADCP, VMP, and CTD data have been processed into relatively small
    % single files (or folders in the case of UHDAS) that span the entire
    % cruise. For these files, load all of them here and assign vessel numbers
    % to each file. Also define a name to use in the output structure for each
    % data source.
    %
    % No data is processed just yet -- we're simply loading existing
    % intermedately processed files/folders and assigning vessel and section
    % numbers for later use. This process can be easily extended to new
    % instrument types or alternative data sources.

    % ----------- UHDAS ADCP data ----------- %
    % Fucent has made some combined .mat files for UHDAS ADCP data. Each ship
    % has a file with variables like this:
    %    wh1200: [1x1 struct]
    %     wh300: [1x1 struct]
    %     wh600: [1x1 struct]
    % wh600_no4: [1x1 struct]

    % We want to load all of these into a single cell array that we can loop
    % through.

    % Locate these files
    uhdas = dir(fullfile(proc_dir,'adcp_ship','*.mat'));

    % Clear any existing loaded adcp data
    clear adcp
    nadcp = 0;

    % Load each file temporarily and insert its structures into the new adcp
    % cell array. Also assign a vessel number and output structure name.
    for i = 1:length(uhdas)
        tmp = load(fullfile(uhdas(i).folder,uhdas(i).name));
        flds = fields(tmp);
        for ii = 1:length(flds) % Looping over wh1200, wh300, etc.
            nadcp = nadcp+1;
            adcp{nadcp} = tmp.(flds{ii}); % insert into combined cell array

            % Assign vessel number and output structure name. These files are
            % named like 'SUNRISE2021_PE_ADCP.mat', so we can identify the
            % vessel with:
            str_parts = strsplit(uhdas(i).name,'_');
            adcp{nadcp}.vnum = find(strcmp(lower({vessels.shortname}),lower(str_parts{2})));
            adcp{nadcp}.name = [vessels(adcp{nadcp}.vnum).shortname '_' flds{ii}];
        end
        clear tmp
    end

    % ----------- VMP data ----------- %
    vmp_files = dir(fullfile(proc_dir,'VMP','SUNRISE2021_*_combo.mat'));
    clear vmp
    for i = 1:length(vmp_files)
        tmp = load(fullfile(vmp_files(i).folder,vmp_files(i).name));
        vmp(i) = tmp.vmp_combo;
        clear tmp;
    end
    % Assign vessel numbers and output structure names to each VMP dataset
    for i = 1:length(vmp)
        % These files are named something like: SUNRISE2021_PE_combo.mat.
        % Extract the short vessel name and find the matching entry in the 'vessels' structure.
        str_parts = strsplit(vmp_files(i).name,'_');
        vmp(i).vnum = find(strcmp(lower({vessels.shortname}), lower(str_parts{2})));
        vmp(i).name = [vessels(vmp(i).vnum).shortname];
    end

    % ----------- CTD data ----------- %
    ctd_files = dir(fullfile(proc_dir,'CTD','SUNRISE2021_*_CTD_combo.mat'));
    clear ctd
    for i = 1:length(ctd_files)
        tmp = load(fullfile(ctd_files(i).folder,ctd_files(i).name));
        ctd(i) = tmp.CTD_combo;
        clear tmp;
    end
    % Assign vessel numbers and output structure names to each CTD dataset
    for i = 1:length(ctd)
        % These files are named something like: SUNRISE2021_PE_combo.mat.
        % Extract the short vessel name and find the matching entry in the 'vessels' structure.
        str_parts = strsplit(ctd_files(i).name,'_');
        ctd(i).vnum = find(strcmp(lower({vessels.shortname}), lower(str_parts{2})));
        ctd(i).name = [vessels(ctd(i).vnum).shortname];
    end

    %% Prepare to load multi-file data sources
    % The T-chain and RHIB ADCP data files are larger and organized by
    % deployment. To avoid loading everything at once, load the time vectors
    % from these files and define for all files:
    % - a vessel number
    % - a file start/end time
    % This will allow us to load data as-needed when looping over sections.
    % Again, this is easily extended to arbitrary deployment-based data
    % structures for instruments/platforms deployed for multiple sections.

    files = struct(); % Contains filepaths and start/end times of data within

    % ----------- RHIB ADCP data ----------- %
    % These are large files, so we'll first determine the time bounds and
    % vessel of each file and then load them as needed.
    files.rhib = dir(fullfile(proc_dir,'adcp_rhib','*','*.mat')); % file metadata
    for i = 1:length(files.rhib)
        % Load time vector from file
        tmp = load(fullfile(files.rhib(i).folder,files.rhib(i).name),'time');

        % Add start and end time to metadata
        files.rhib(i).start_time = min(tmp.time);
        files.rhib(i).end_time = max(tmp.time);

        % Figure out which vessel each file belongs to and store the vessel index
        [~,vname] = fileparts(files.rhib(i).folder);
        files.rhib(i).vnum = find(strcmp(lower({vessels.name}),lower(vname)));
    end

    % ----------- TChain data ----------- %
    % Similar to the RHIB files.
    files.tchain = dir(fullfile(proc_dir,'tchain','*.mat')); % file metadata
    for i = 1:length(files.tchain)
        % Load time vector from file
        tmp = load(fullfile(files.tchain(i).folder,files.tchain(i).name),'dn');

        % Add start and end time to metadata
        files.tchain(i).start_time = min(tmp.dn);
        files.tchain(i).end_time = max(tmp.dn);

        % Figure out which vessel each file belongs to and store vessel index
        files.tchain(i).vnum = find(cellfun(@(s) contains(files.tchain(i).name,s), {vessels.name}));
    end

    % ----------- Done gathering multi-section file metadata -----------%

    % Use the multi-file data source information defined above to assign file
    % numbers to each vessel & section

    % Create section start/end time matrix
    % NaN's if vessel not in section
    sec_start = nan(max(secs.n),length(vessels));
    sec_end   = nan(max(secs.n),length(vessels));
    for v = 1:length(vessels)
        % Filter sections to
        vsecs = secs(find(strcmp(secs.vessel,vessels(v).name)),:);
        for i = 1:height(vsecs)
            sec_start(vsecs.n(i),v) = datenum(vsecs.start_utc(i));
            sec_end(vsecs.n(i),v) = datenum(vsecs.end_utc(i));
        end
    end

    % Map instrument file numbers to sections for all vessels
    insts = {'tchain','rhib'};
    inst_fullnames = {'TChain','RHIB ADCP'};
    for ii = 1:length(insts)
        nfile.(insts{ii}) = zeros(max(secs.n),length(vessels));
        for v = 1:length(vessels) % iterate over vessels
            vessels(v).file_loaded.(insts{ii}) = 0;  % initialize currently loaded filenum
            vessels(v).data_loaded.(insts{ii}) = []; % initialize currently loaded data
            for i = 1:max(secs.n)
                % Find a file that:
                % 1) belongs to the vessel
                % 2) starts before the section start time
                % 3) ends after the section end time
                parent_fnum = find([files.(insts{ii}).vnum]       == v & ...
                                   [files.(insts{ii}).start_time] <= sec_start(i,v) & ...
                                   [files.(insts{ii}).end_time]   >= sec_end(i,v));
                if ~isempty(parent_fnum)
                    % If such a file is found, store it
                    nfile.(insts{ii})(i,v) = parent_fnum;
                end
            end
        end
    end

    %% Generate section files for sections in survey
    % This is the main processing loop. We use all of the metadata above to load
    % files as needed, index by time for every section, and save an output file.

    %========================================================================
    % The sections handling each instrument type can be modified as needed to
    % satisfy any desired naming conventions. I simply took existing
    % intermediate data and indexed by time while preserving variable names.
    %========================================================================

    for i = 1:max(secs.n)
        section = struct();

        fprintf('├── Section %d\n',i);

        %% UHDAS ADCP data
        for a = 1:length(adcp)

            % FIXME: All processed files should have a datenum field. Workaround
            % for now.
            if ~isfield(adcp{a},'dn');
                adcp{a}.dn = datenum(adcp{a}.time');
            end

            % Check if this ADCP contains any data during the section
            idx = adcp{a}.dn >= sec_start(i,adcp{a}.vnum) & ...
                  adcp{a}.dn <= sec_end(i,adcp{a}.vnum);

            % If it does, extract the indices contained in the section.
            if sum(idx)>0
                fprintf('   ├── UHDAS ADCP in section: %s\n',adcp{a}.name);
                flds = setdiff(fields(adcp{a}),{'vnum','name'});
                for f = 1:length(flds)
                    % Variables with a time dimension get indexed and stored
                    if size(adcp{a}.(flds{f}),2) == length(idx)
                        section.adcp.(adcp{a}.name).(flds{f}) = adcp{a}.(flds{f})(:,idx);
                    else
                        % Variables with no time dimension (e.g. depth) get stored in full
                        section.adcp.(adcp{a}.name).(flds{f}) = adcp{a}.(flds{f});
                    end
                end
            end
        end

        %% VMP data
        for v = 1:length(vmp)
            % Check if this VMP structure contains any data during the section
            idx = vmp(v).dn >= sec_start(i,vmp(v).vnum) & ...
                  vmp(v).dn <= sec_end(i,vmp(v).vnum);

            % If it does, extract the indices contained in the section.
            if sum(idx)>0
                fprintf('   ├── VMP in section: %s\n',vmp(v).name);
                flds = setdiff(fields(vmp(v)),{'name','vnum'});
                for f = 1:length(flds)
                    % Variables with a time dimension get indexed and stored
                    if size(vmp(v).(flds{f}),2) == length(idx)
                        section.vmp.(vmp(v).name).(flds{f}) = vmp(v).(flds{f})(:,idx);
                    else
                        % Variables with no time dimension (e.g. depth) get stored
                        section.vmp.(vmp(v).name).(flds{f}) = vmp(v).(flds{f});
                    end
                end
            end
        end

        %% CTD data
        for v = 1:length(ctd)
            % Check if this CTD structure contains any data during the section
            idx = ctd(v).dn >= sec_start(i,ctd(v).vnum) & ...
                  ctd(v).dn <= sec_end(i,ctd(v).vnum);

            % If it does, extract the indices contained in the section.
            if sum(idx)>0
                fprintf('   ├── CTD in section: %s\n',ctd(v).name);
                flds = setdiff(fields(ctd(v)),{'name','vnum'});
                for f = 1:length(flds)
                    % Variables with a time dimension get indexed and stored
                    if size(ctd(v).(flds{f}),2) == length(idx)
                        section.ctd.(ctd(v).name).(flds{f}) = ctd(v).(flds{f})(:,idx);
                    else
                        % Variables with no time dimension (e.g. depth) get stored
                        section.ctd.(ctd(v).name).(flds{f}) = ctd(v).(flds{f});
                    end
                end
            end
        end


        %% Intermediate data organized by deployments (TChains, RHIB ADCP)
        for ii = 1:length(insts)
            for v = 1:length(vessels)

                % ------------ Load data as needed ------------ %

                % Find a file containing data from this section
                fnum = nfile.(insts{ii})(i,v);

                % If one exists, convert its data to the desired format
                if fnum > 0

                    % Report some information
                    name = vessels(v).shortname; % name of structure in output
                    fprintf('   ├── %s in section: %s\n',inst_fullnames{ii},name);

                    % Load new data if necessary
                    if vessels(v).file_loaded.(insts{ii}) ~= fnum;
                        fprintf('     ├── Loading %s...',files.(insts{ii})(fnum).name);
                        fname = fullfile(files.(insts{ii})(fnum).folder,files.(insts{ii})(fnum).name);
                        vessels(v).data_loaded.(insts{ii}) = load(fname);
                        vessels(v).file_loaded.(insts{ii}) = fnum;
                        fprintf('\r     ├── Loaded %s\n',files.(insts{ii})(fnum).name);
                    end

                    %% Index by time and extract/rename desired data from intermediate files
                    dat =@(inst) vessels(v).data_loaded.(inst); % Access loaded instrument file
                    switch insts{ii}

                      case 'rhib'
                        % Identify time indices in section
                        idx = dat(insts{ii}).time >= sec_start(i,v) & ...
                              dat(insts{ii}).time <= sec_end(i,v);
                        if sum(idx) > 0
                            % Extract/rename data
                            section.adcp.(name).dn    = dat(insts{ii}).nuc_time(idx);
                            section.adcp.(name).depth = dat(insts{ii}).cell_depth;
                            section.adcp.(name).z     = -dat(insts{ii}).cell_depth;
                            section.adcp.(name).u     = squeeze(dat(insts{ii}).vel(:,1,idx));
                            section.adcp.(name).v     = squeeze(dat(insts{ii}).vel(:,2,idx));
                            section.adcp.(name).w     = squeeze(dat(insts{ii}).vel(:,3,idx));
                        end

                      case 'tchain'
                        % Identify time indices in section
                        idx = dat(insts{ii}).dn >= sec_start(i,v) & ...
                              dat(insts{ii}).dn <= sec_end(i,v);
                        if sum(idx) > 0;
                            % Extract/rename data
                            section.tchain.(name).dn  = dat(insts{ii}).dn(idx);
                            section.tchain.(name).lat = dat(insts{ii}).lat(idx);
                            section.tchain.(name).lon = dat(insts{ii}).lon(idx);
                            section.tchain.(name).z   = dat(insts{ii}).z(:,idx);
                            section.tchain.(name).t   = dat(insts{ii}).t(:,idx);
                            section.tchain.(name).p   = dat(insts{ii}).p(:,idx);
                            section.tchain.(name).s   = dat(insts{ii}).s(:,idx);
                            section.tchain.(name).pos = dat(insts{ii}).pos;
                        end
                    end % convert data

                end % if data exists
            end % loop over vessels
        end % loop over instruments

        % Load sensor offsets
        % Compute position of all instruments
        % Add to data structures

        % Save output
        file_name = sprintf('SUNRISE_2021_survey_%02d_section_%02d.mat',s,i);
        save(fullfile(dir_out,file_name),'-struct','section');

        disp(['Saved ' file_name]);

    end % loop over sections
end
