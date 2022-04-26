%% initial setup
clear
close all

%%% setting path and pre-loading files
addpath('../_Config')
Process_Mode = 'Biosonics';
data_path %% all data path and library

%%
DT4_list = dir(Bioson_RAWD_Path);
DT4_list = DT4_list(endsWith({DT4_list.name},'.dt4'));

for i = 1:length(DT4_list)
    pings = readbio([DT4_list(i).folder '/' DT4_list(i).name],1,10,5);
    save([Bioson_RAWM_Path '/' DT4_list(i).name(1:end-4) '.mat'],'-struct','pings')
    
    movefile([DT4_list(i).folder '/' DT4_list(i).name],[DT4_list(i).folder '/done/' DT4_list(i).name])
end