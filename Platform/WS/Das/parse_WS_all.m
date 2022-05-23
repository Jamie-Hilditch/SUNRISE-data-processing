clear all, close all

% List all WS DAS folders
folders = dir('../../data/raw/ship_DAS/WaltonSmith/');
[~,idx] = setdiff({folders.name},{'.','..'});
folders = folders(idx);

file_out = '../../data/intermediate/ship_DAS/WS_DAS.mat';

% In most cases, read all files. In some cases, read specific files:
match = {'Depth'            , 'SDDBT';
         'PortRMYoung'      , 'PTU';
         'RADD1'            , 'LN2';
         'RADD2'            , 'LN1';
         'SpeedLog'         , 'Water';
         'StarboardRMYoung' , 'Wind';
         'Turner C7'        , 'C3'};

% Skip these directories. Mostly extra GPS output with different filetypes. Will
% need to parse these separately (see parse_WS_gps.m) if we're interested in
% this data.
skip = {'GPS1','GPS2','GyroCompass','Lister','POSMV'};

out = struct();
for i = 1:length(folders)

    % Create name for output structure
    name = strrep(folders(i).name,' ','_');
    name = strrep(name,'-','_');
    if ~ismember(folders(i).name,skip)
        fprintf('Parsing %s to %s\n',folders(i).name,name);
        dir_in = fullfile(folders(i).folder,folders(i).name);

        if ismember(folders(i).name,match(:,1));
            pat = match{find(strcmp(folders(i).name,match(:,1))),2};
        else
            pat = '';
        end

        out = struct();
        out.(name) = parse_WS_DAS(dir_in,pat);
        file_out = sprintf('../../data/intermediate/ship_DAS/WS_%s.mat',name);
        save(file_out,'-v7.3','-struct','out');
        disp(['Saved ', file_out])

    end
end
