clear
close

Processed_Path = '../../../../data/Processed/';

process_mat_path = dir(fullfile(Processed_Path, '**/*.*'));
Processed_Path = [process_mat_path(1).folder '/'];
process_mat_path = process_mat_path(startsWith({process_mat_path.name},'SUNRISE2021_') & endsWith({process_mat_path.name},'_Processed.mat'));

for i = 1:length(process_mat_path)
    nc_folder = [Processed_Path(1:end-1) '_NC/' process_mat_path(i).folder(37:end)];
    nc_path = [nc_folder '/' process_mat_path(i).name(1:end-3) 'nc'];
    mat_path = [process_mat_path(i).folder '/' process_mat_path(i).name];
    if ~exist(nc_folder,'dir')
        mkdir(nc_folder)
    end
    
    
    var_list = whos('-file',mat_path);
    var_list = var_list(strcmp({var_list.class},'double'));
    
    dn_idx = find(strcmp({var_list.name},'dn'));
    dn_length = max(var_list(dn_idx).size);
    
    depth_idx = find(strcmp({var_list.name},'depth'));
    try
        depth_length = max(var_list(depth_idx).size);
        if depth_length == dn_length
            depth_length = min(var_list(depth_idx).size);
        end
    catch
        depth_length = 0;
        disp(['No depth data:' process_mat_path(i).name(1:end-3) 'nc'])
    end
    
    
    mat_data = matfile(mat_path);
    if exist(nc_path,'file')
        %% old nc append
        
        old_dn_length = ncinfo(nc_path,'dn');
        old_dn_length = old_dn_length.Dimensions.Length;
        
        if dn_length > old_dn_length
            for j = 1:length(var_list)
                var_size = var_list(j).size;
                if ndims(var_size)<=2
                    dn_dim = find(var_size == dn_length);
                    depth_dim = find(var_size == depth_length);
                    if isempty(depth_dim)
                        other_length = var_size;
                        other_length(dn_dim) = [];
                    else
                        other_length = [];
                    end
                    
                    if ~isempty(dn_dim.*depth_dim)
                        if dn_dim == 2
                            ncwrite(nc_path,var_list(j).name,mat_data.(var_list(j).name),[1 old_dn_length+1])
                        else                            
                            ncwrite(nc_path,var_list(j).name,mat_data.(var_list(j).name),[old_dn_length+1 1])                            
                        end
                    elseif ~isempty(dn_dim)
                        if other_length==1
                            ncwrite(nc_path,var_list(j).name,mat_data.(var_list(j).name),old_dn_length+1)
                        else
                            if dn_dim == 2
                                ncwrite(nc_path,var_list(j).name,mat_data.(var_list(j).name),[1 old_dn_length+1])
                            else
                                ncwrite(nc_path,var_list(j).name,mat_data.(var_list(j).name),[old_dn_length+1 1])
                            end
                        end
                    end
                    
                else
                    disp(['Skip ' var_list(j).name '(' process_mat_path(i).name(1:end-3) 'nc)'])
                end
            end
        else
            disp(['No New Data: '  process_mat_path(i).name(1:end-3) 'nc'])
        end
    else
        %% new file create
        for j = 1:length(var_list)
            var_size = var_list(j).size;
            if ndims(var_size)<=2
                dn_dim = find(var_size == dn_length);
                depth_dim = find(var_size == depth_length);
                if isempty(depth_dim)
                    other_length = var_size;
                    other_length(dn_dim) = [];
                    
                else
                    other_length = [];
                end
                
                if ~isempty(dn_dim.*depth_dim)
                    if dn_dim == 2
                        nccreate(nc_path,var_list(j).name, ...
                            'Dimensions',{'depth',depth_length,'dn',inf},...
                            'Format','netcdf4')
                        ncwrite(nc_path,var_list(j).name,mat_data.(var_list(j).name))
                    else
                        nccreate(nc_path,var_list(j).name, ...
                            'Dimensions',{'dn',inf,'depth',depth_length},...
                            'Format','netcdf4')
                    end
                    
                elseif ~isempty(dn_dim)
                    if other_length==1
                        nccreate(nc_path,var_list(j).name, ...
                            'Dimensions',{'dn',inf},...
                            'Format','netcdf4')
                        ncwrite(nc_path,var_list(j).name,mat_data.(var_list(j).name))
                    else
                        if dn_dim == 2
                            nccreate(nc_path,var_list(j).name, ...
                                'Dimensions',{['other_' var_list(j).name],other_length,'dn',inf},...
                                'Format','netcdf4')
                            ncwrite(nc_path,var_list(j).name,mat_data.(var_list(j).name))
                        else
                            nccreate(nc_path,var_list(j).name, ...
                                'Dimensions',{'dn',inf,['other_' var_list(j).name],other_length},...
                                'Format','netcdf4')
                            ncwrite(nc_path,var_list(j).name,mat_data.(var_list(j).name))
                        end
                    end
                end
                
            else
                disp(['Skip ' var_list(j).name '(' process_mat_path(i).name(1:end-3) 'nc)'])
            end
        end
    end
end



