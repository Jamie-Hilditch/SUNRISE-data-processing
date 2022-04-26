%% initial setup
clear
close all

%%% setting path and pre-loading files
addpath('../_Config')
Process_Mode = 'Biosonics';
data_path %% all data path and library

MAT_list = dir(Bioson_RAWM_Path);
MAT_list = MAT_list(endsWith({MAT_list.name},'.mat'));

BS_name = {'datenum','lat','lon','sample','Sv','TS'};
BS_new_name = {'dn','lat','lon','BS_sample','BS_Sv','BS_TS'};

for i = 1:length(MAT_list)
    
    BS_data = load([MAT_list(i).folder '/' MAT_list(i).name]);
    
    
    if exist('BS_combo_temp','var')
        N = length(BS_combo_temp.dn);
    end
    N_chunk = length(BS_data.datenum);
    
    for k = 1:length(BS_name)
        temp = BS_data.(BS_name{k});
        if i >1
            
            if isvector(temp)
                BS_combo_temp.(BS_new_name{k})(N+1:N+N_chunk) = temp;
            else
                BS_combo_temp.(BS_new_name{k})(:,N+1:N+N_chunk) = temp;
            end
        else
            BS_combo_temp.(BS_new_name{k}) = temp;
        end
    end
    movefile([MAT_list(i).folder '/' MAT_list(i).name],[MAT_list(i).folder '/done/' MAT_list(i).name])
end
BS_combo_temp.depth = BS_data.depth;


if exist([Bioson_PROC_final_Path Prefix '_Biosonics_Processed.mat'],'file')
    BS_combo = load([Bioson_PROC_final_Path Prefix '_Biosonics_Processed.mat']);
    N = length(BS_combo.dn);
    N_chunk = length(BS_combo_temp.dn);
    
    for k = 1:length(BS_new_name)
        temp = BS_combo_temp.(BS_new_name{k});
        if isvector(temp)
            BS_combo.(BS_new_name{k})(1,N+1:N+N_chunk) = temp;
        else
            BS_combo.(BS_new_name{k})(:,N+1:N+N_chunk) = temp;
        end
    end
    
    [~,idx_temp] = sort(BS_combo.dn);
    duplicate_idx = find(diff(BS_combo.dn(1,idx_temp))==0);
    
    if ~isempty(duplicate_idx)
        idx_temp(duplicate_idx) = [];
        for k = 1:length(BS_new_name)
            if isvector(temp)
                BS_combo.(BS_new_name{k}) = BS_combo.(BS_new_name{k})(1,idx_temp);
            else
                BS_combo.(BS_new_name{k}) = BS_combo.(BS_new_name{k})(:,idx_temp);
            end
        end
    end
    save([Bioson_PROC_final_Path Prefix '_Biosonics_Processed.mat'],'-v7.3','-struct','BS_combo')
        
else
    BS_combo = BS_combo_temp;
    save([Bioson_PROC_final_Path Prefix '_Biosonics_Processed.mat'],'-v7.3','-struct','BS_combo')
end

