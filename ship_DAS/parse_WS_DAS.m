%% parse_WS_DAS.m
% Usage: parse_WS_DAS(dir_in, name)
% Description: Parse data from any files in the Walton Smith's shipboard DAS folder.
%              This will create a table containing data from all .dat files of the same type.
% Inputs: dir_in   - path to the directory containing files
%         name     - File type to parse, i.e. all files named "*[name].dat"
%
% Outputs: tbl - Table containing data from all files of the named type
%
% Author: Dylan Winters
% Created: 2022-01-06

function tbl = parse_WS_DAS(dir_in,name)

if isempty(name)
    files = dir(fullfile(dir_in,['*.dat']));
else
    files = dir(fullfile(dir_in,['*' name '.dat']));
end

% Create import options; explicitly define date format
fname = fullfile(files(1).folder,files(1).name);
opts = detectImportOptions(fname);
if ismember('ComputerDate',opts.VariableNames)
    opts = setvaropts(opts,'ComputerDate','InputFormat','MM/dd/uuuu');
elseif strcmp(opts.VariableTypes{1},'datetime')
    opts = setvaropts(opts,opts.VariableNames{1},'InputFormat','MM/dd/uuuu');
end


% MATLAB will warn when it has to rename columns (e.g. removing spaces).
% Disable this.
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

% Read all files
tbl = [];
for i = 1:length(files)
    fname = fullfile(files(i).folder,files(i).name);

    % Read first file, concat all subsequent files
    if isempty(tbl)
        tbl = readtable(fname,opts);
    else
        try
            tbl = cat(1,tbl,readtable(fname,opts));
        catch err
            fprintf('\nSkipped %s\n',fname)
        end
    end
    fprintf('\rRead %s [%d of %d]',files(i).name,i,length(files));
end
fprintf('\n')
