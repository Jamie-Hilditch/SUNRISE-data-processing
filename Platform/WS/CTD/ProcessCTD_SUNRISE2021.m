%%% This Code is rewritten by Fucent in 2021 to process VMP data at the field - 2021 SUNRISE.
%%% It does all the processing procedure to generate dissipation rate
%%% It based on quick_look and OSU 2019 version.
%%% code start %%%
%% initial setup
clear
close all

%%% setting path and pre-loading files
addpath('../_Config')
Process_Mode = 'CTD';
data_path %% all data path and library

%% version check
% if odas_version_info ~= 4.4
%     error('ERROR: please use odas v4.4')
% end

%% vmp log
logid = fopen('./ctd_process_log.txt','a');


%% reading file
Raw_list = dir([CTD_RAWR_Path '*.rsk']);
%Raw_list = Raw_list(~endsWith({Raw_list.name}, '_original.P'));

proc_idx = 1:2;%length(Raw_list);

for i = proc_idx
    if exist([CTD_RAWM_Path Raw_list(i).name(1:end-4) '.mat'],'file')
        raw_mat = load([CTD_RAWM_Path Raw_list(i).name(1:end-4) '.mat']);
    else
        raw_mat = RSKopen([CTD_RAWR_Path Raw_list(i).name]);
        raw_mat = RSKreaddata(raw_mat);
        
        [raw_mat,~] = RSKcorrecthold(raw_mat);
        %RSKdespike %% should do it after see actual profile
        %RSKsmooth %% should do it after see actual profile
        [~,temp1] = min(abs(raw_mat.data.tstamp - raw_mat.data.tstamp(1) - 1/1440));
        raw_mat = RSKderiveseapressure(raw_mat,'patm',nanmean(raw_mat.data.values(1:temp1,getchannelindex(raw_mat, 'Pressure'))));
        raw_mat = RSKderivedepth(raw_mat,'latitude',lat_c);
        raw_mat = RSKderivevelocity(raw_mat);
        
        CT_lag = RSKcalculateCTlag(raw_mat);
        raw_mat = RSKalignchannel(raw_mat, 'channel', 'Conductivity', 'lag', CT_lag);
        
        raw_mat = RSKderivesalinity(raw_mat);
        raw_mat = RSKderivetheta(raw_mat,'latitude',lat_c,'longitude',lon_c);
        raw_mat = RSKderiveSA(raw_mat,'latitude',lat_c,'longitude',lon_c);
        raw_mat = RSKderivesigma(raw_mat,'latitude',lat_c,'longitude',lon_c);
        if isempty(find(ismember({raw_mat.channels.longName},'Dissolved O2'), 1))
            temp.values = raw_mat.data.tstamp*nan;
            
            raw_mat = RSKaddchannel(raw_mat,'data',temp,'channel','Dissolved O2','unit','%');
            raw_mat = RSKaddchannel(raw_mat,'data',temp,'channel','Dissolved O22','unit','ml/l');
        else
            raw_mat = RSKderiveO2(raw_mat, 'toDerive','concentration', 'unit', 'ml/l');
        end
        
        save([CTD_RAWM_Path Raw_list(i).name(1:end-4) '.mat'],'-struct','raw_mat','-v7.3')
        movefile([CTD_RAWR_Path Raw_list(i).name],[CTD_RAWR_Path 'done/' Raw_list(i).name])
        fprintf(logid,'%19s %8s Converting RSK-file.... \n',datestr(now,'yyyy/mm/dd HH:MM:SS'),Raw_list(i).name);
    end
    
    profile = RSKtimeseries2profiles(raw_mat,'pressureThreshold',3);
    profile = RSKremoveloops(profile);    
    profile = RSKtrim(profile, 'reference', 'sea pressure', 'range',[-1 0.5], 'action', 'remove');    
    
    profile_idx = find(ismember({profile.data.direction},'down'));            
    %%
    for j = 1:length(profile_idx)
        CTD_name = {'Pressure','Depth','Conductivity','Salinity','Absolute Salinity','Temperature','Potential Temperature',...
            'Density Anomaly','Dissolved O2','Dissolved O22'};
        CTD_short = {'P','depth','C','SP','SA','T','theta','sigma','DO2R','DO2A'};
        
        v_idx = find(~isnan(profile.data(profile_idx(j)).values(:,3)));
        N = length(v_idx);
        
        for k = 1 :length(CTD_name)
            idx = find(ismember({profile.channels.longName},CTD_name{k}));
            if j == 1
                eval(['ctd.' CTD_short{k} '= nan(1000,length(profile_idx));']);
                eval(['ctd.' CTD_short{k} '(1:N,j) = profile.data(profile_idx(j)).values(v_idx,idx);']);
            elseif j == length(profile_idx)
                eval(['ctd.' CTD_short{k} '(1:N,j) = profile.data(profile_idx(j)).values(v_idx,idx);']);
                if k==1
                    eval(['ctd.' CTD_short{k} '(find(isnan(nanmean(ctd.P,2)),1):end,:) = [];']);
                else
                    eval(['ctd.' CTD_short{k} '(size(ctd.P,1)+1:end,:) = [];']);
                end
            else
                eval(['ctd.' CTD_short{k} '(1:N,j) = profile.data(profile_idx(j)).values(v_idx,idx);']);
            end
        end
        
        ctd.dn(j) = profile.data(profile_idx(j)).tstamp(v_idx(1));
        
        [~,time_idx] = min(abs(ship.dn-ctd.dn(j)));
        ctd.lat(j) = ship.lat(time_idx,1);
        ctd.lon(j) = ship.lon(time_idx,1);
        ctd.dist_ctd(j) = ship.dist_ship(time_idx,1);
        ctd.data_num(j) = profile.epochs.startTime+profile.instruments.serialID*1e6; % RBR Serial ID
        ctd.SN(j) = profile.instruments.serialID;
        ctd.profile(j) = j;
    end
    save([CTD_PROC_Path Prefix '_raw_ctd' datestr(profile.epochs.startTime,'_yyyymmddHHMMSS') '.mat'],'-struct','ctd','-v7.3')
    clear ctd
end
fclose(logid);
