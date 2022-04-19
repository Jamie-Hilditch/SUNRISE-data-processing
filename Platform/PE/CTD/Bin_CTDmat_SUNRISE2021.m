%%% This Code is rewritten by Fucent in 2021 to process VMP data at the field - 2021 SUNRISE.
%%% Raw combine VMP data is binned in 1 m in this code.
%%% generate matfile and NETCDF for the transfer purpose.
%%% code start %%%

%% initial setup
clear
close all

%%%
level = 40;

%%% setting path and pre-loading files
addpath('../../general toolbox/')
dictionary = 'ctd';
data_path %% all data path and library

%% reading file
Comb_mat_list = dir([CTD_PROC_final_Path '/' basename '_raw_ctd*.mat']);
Comb_timelist = cell2mat(extractfield(Comb_mat_list,'name')');
Comb_timelist = str2double(cellstr(char(Comb_timelist(:,16:end-4))));
[Comb_timelist,comb_file_idx] = sort(Comb_timelist);
Comb_mat_list = Comb_mat_list(comb_file_idx);

%% what will be binned
ctd_names = {'C','SP','SA','T','theta','sigma','O2R','O2A'};

for file_i = 1:length(Comb_mat_list)
    ctd = load([Comb_mat_list(file_i).folder '/' Comb_mat_list(file_i).name]);
    
%% create matrix
for i = 1:length(ctd_names)
    eval(['temp = ctd.' char(ctd_names(i)) ';']);
    
    temp3 = nan(level,size(ctd.depth,2));
    for z = 1:level
        mask = nan(size(ctd.depth));
        mask(ctd.depth >= (z-1) & ctd.depth < z) = 1;
        temp3(z,:) = nanmean(mask.*temp);
    end
    temp3(temp3==0) = nan;
    
    eval(['CTD_combo_temp.' char(ctd_names(i)) ' = temp3;'])
end

CTD_combo_temp.dn = ctd.dn;
CTD_combo_temp.lat = ctd.lat;
CTD_combo_temp.lon = ctd.lon;
CTD_combo_temp.dist = ctd.dist;
CTD_combo_temp.data_num = ctd.data_num;
CTD_combo_temp.profile = ctd.profile;
CTD_combo_temp.depth = (-0.5:-1:-level+0.5);


if exist([CTD_PROC_final_Path '/' basename '_CTD_combo.mat'],'file')
    load([CTD_PROC_final_Path '/' basename '_CTD_combo.mat'])
    
    %oldtime_idx = find(CTD_combo.dn == CTD_combo_temp.dn(1));
    [~,oldtime_idx] = min(abs(CTD_combo.dn - CTD_combo_temp.dn(1)));
    
    if oldtime_idx == length(CTD_combo.dn)
        N_str = length(CTD_combo.dn)+1;
    else
        N_str = oldtime_idx+1;
    end
    
    ctd_names = {'C','SP','SA','T','theta','sigma','O2R','O2A'};
    
    for i = 1:length(ctd_names)
        eval(['CTD_combo.' char(ctd_names(i)) '(:,N_str:N_str-1+length(CTD_combo_temp.dn)) = CTD_combo_temp.' char(ctd_names(i)) ';'])
    end
    
    CTD_combo.dn(N_str:N_str-1+length(CTD_combo_temp.dn)) = CTD_combo_temp.dn;
    CTD_combo.lat(N_str:N_str-1+length(CTD_combo_temp.dn)) = CTD_combo_temp.lat;
    CTD_combo.lon(N_str:N_str-1+length(CTD_combo_temp.dn)) = CTD_combo_temp.lon;
    CTD_combo.dist(N_str:N_str-1+length(CTD_combo_temp.dn)) = CTD_combo_temp.dist;
    CTD_combo.data_num(N_str:N_str-1+length(CTD_combo_temp.dn)) = CTD_combo_temp.data_num;
    CTD_combo.profile(N_str:N_str-1+length(CTD_combo_temp.dn)) = CTD_combo_temp.profile;
    CTD_combo.depth = (-0.5:-1:-level+0.5);
    
else
    CTD_combo = CTD_combo_temp;
end

save([CTD_PROC_final_Path '/' basename '_CTD_combo.mat'],'CTD_combo','-v7.3')
end




