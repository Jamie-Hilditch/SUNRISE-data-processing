%% initial setup
clear
close all

%%% setting path and pre-loading files
addpath('../_Config')
Process_Mode = 'Tchain';
data_path %% all data path and library

%%%

dir_raw = dir(TCn_DATA_Path);
dir_raw = dir_raw(startsWith({dir_raw.name},'deploy'));


warning off
for i = 1:length(dir_raw)
    if exist([dir_raw(i).folder '/' dir_raw(i).name '/instrument_depths.csv'],'file')
        if exist([dir_raw(i).folder '/' dir_raw(i).name '/instrument_depths_og.csv'],'file')
            t = readtable([dir_raw(i).folder '/' dir_raw(i).name '/instrument_depths_og.csv']);
            [t.(dir_raw(i).name),t.SectionDuration,t.zero_pressure_interval,t.Skip] = deal(nan(size(t.depth_m_)));
            [t.SectionDuration] = deal(NaT(size(t.depth_m_)));
            t = movevars(t,6:8,'Before',1);
            t = movevars(t,9,'Before',8);
            
            t_skip = readtable([dir_raw(i).folder '/' dir_raw(i).name '/instrument_depths.csv']);
            
            idx = find(~ismember(t.serialnum,t_skip.serialnum));
            t.Skip(idx) = 1;
            
        else
            t = readtable([dir_raw(i).folder '/' dir_raw(i).name '/instrument_depths.csv']);
            [t.(dir_raw(i).name),t.SectionDuration,t.zero_pressure_interval,t.Skip] = deal(nan(size(t.depth_m_)));
            [t.SectionDuration] = deal(NaT(size(t.depth_m_)));
            t = movevars(t,6:8,'Before',1);
            t = movevars(t,9,'Before',8);
        end
        
        t_time = readtable([dir_raw(i).folder '/' dir_raw(i).name '/section_start_end_time.csv']);
        t.SectionDuration(1:2) = [t_time.start_time t_time.end_time];
        t.var = nan(size(t.depth_m_));
        
        for j = 2:width(t)
            t.Properties.VariableNames{j} = [t.Properties.VariableNames{j} int2str(i)];
        end
        
        if i ==1
            if exist('./Deployment_Info.cvs','file')
                t_final = readtable('./Deployment_Info.cvs');
            else
                t_final = t;
            end
        else
            if height(t_final)>height(t)
                
                t(height(t)+1:height(t_final),1) = {nan};
                
                t_final(:,end+1:end+10) = t;
                
            elseif  height(t_final)<height(t)
                t_final(height(t_final)+1:height(t),1) =  {nan};
                t_final(:,end+1:end+10) = t;
                
            else
                t_final(:,end+1:end+10) = t;
            end
            t_final.Properties.VariableNames(end-9:end) = t.Properties.VariableNames;
            % t_final = outerjoin(t_final,t);
        end
    end
end
warning on

writetable(t_final,'Deployment_Info.csv')
%t_output =