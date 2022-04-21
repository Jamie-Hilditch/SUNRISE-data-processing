%%% This Code is rewritten by Fucent in 2021 to process VMP data at the field - 2021 SUNRISE.
%%% It does all the processing procedure to generate dissipation rate
%%% It based on quick_look and OSU 2019 version.
%%% code start %%%
%% initial setup
clear
close all

%%% setting path and pre-loading files
addpath('../_Config')
Process_Mode = 'VMP';
data_path %% all data path and library

%% version check
if odas_version_info ~= 4.4
    error('ERROR: please use odas v4.4')
end

%% vmp log
logid = fopen('./vmp_process_log.txt','a');

%% parameter
%%% despike
% [thresh, smooth, and length] (in seconds) -> Rockland default value,
% subject to be altered.

% for shear probes
despike_sh.thres  = 8;
despike_sh.smooth = 0.5;
despike_sh.N_FS = 0.05;

% for piezo-accelerometers
despike_A.thres  = 8;
despike_A.smooth = 0.5;
despike_A.N_FS = 0.05;

%%% dissipation
diss_info.fft_length_sec    = 0.5; % unit second, falling rate~1.2m/s
diss_info.diss_length_fac   = 2; % no unit, >2*fft_length Rec:>3

%% module
bot_mod = 1;

%% bin
level = 40; 
dz = 1; % m

%% reading file
Raw_list = dir([VMP_RAWP_Path '*.P']);
Raw_list = Raw_list(~endsWith({Raw_list.name}, '_original.P') & ~startsWith({Raw_list.name}, '._'));

proc_idx = 1:length(Raw_list);


for i = proc_idx
    if exist([VMP_RAWM_Path Raw_list(i).name(1:end-2) '.mat'],'file')
        raw_mat = load([VMP_RAWM_Path Raw_list(i).name(1:end-2) '.mat']);
    else
        raw_mat = odas_p2mat([VMP_RAWP_Path Raw_list(i).name]);
        save([VMP_RAWM_Path Raw_list(i).name(1:end-2) '.mat'],'-struct','raw_mat','-v7.3')
        fprintf(logid,'%19s %8s Converting P-file.... \n',datestr(now,'yyyy mm/dd HH:MM:SS'),Raw_list(i).name);
    end
    
    %% seperate profile
    profile_idx = get_profile(raw_mat.P_slow, raw_mat.W_slow, 0.5, 0.3, 'down', 10, raw_mat.fs_slow);
    if isempty(profile_idx)
        movefile([VMP_RAWP_Path Raw_list(i).name],[VMP_RAWP_Path '/junk/' Raw_list(i).name])
        fprintf(logid,'%19s %8s No Profile - moving to junk \n',datestr(now,'yyyy/mm/dd HH:MM:SS'),Raw_list(i).name);
    else
        %% checking FP07 thermistor calibration
        cali = 1;
        
        while cali
            %T1
            [T1_0, beta_T1, Lag_T1] = cal_FP07_acc(raw_mat,profile_idx, 'JAC_T','T1',1);
            %T1 coef in cfg file
            T1_0_def = str2double(setupstr(raw_mat.cfgobj, 'T1' , 't_0'));
            B1_0_def = str2double(setupstr(raw_mat.cfgobj, 'T1' , 'beta_1'));
            
            %T2
            [T2_0, beta_T2, Lag_T2] = cal_FP07_acc(raw_mat,profile_idx, 'JAC_T','T2',1);
            %T2 coef in cfg file
            T2_0_def = str2double(setupstr(raw_mat.cfgobj, 'T2' , 't_0'));
            B2_0_def = str2double(setupstr(raw_mat.cfgobj, 'T2' , 'beta_1'));
            
            if abs(T1_0-T1_0_def)>1e-3 || abs(B1_0_def- beta_T1(1)>1e-2) || abs(T2_0-T2_0_def)>1e-3 || abs(B2_0_def- beta_T2(1)>1e-2)
                temp = raw_mat.setupfilestr;
                if abs(T1_0-T1_0_def)>1e-3 || abs(B1_0_def- beta_T1(1)>1e-2)
                    temp = regexprep(temp,sprintf('(SN(.){0,11}= %5s(((.|\\n){0,100})%7.3f))',char(setupstr(raw_mat.cfgobj, 'T1','SN')),str2double(setupstr(raw_mat.cfgobj, 'T1','T_0'))),...
                        sprintf(['SN          = %5s\n'...
                        'beta_1      = %7.2f\n'...
                        'T_0         = %7.3f\n'],...
                        char(setupstr(raw_mat.cfgobj, 'T1','SN')),beta_T1(1),T1_0));
                    fprintf(logid,'%19s %8s T1 Calibration T_0: %7.3f  Beta1: %7.2f   \n',datestr(now,'yyyy/mm/dd HH:MM:SS'),Raw_list(i).name,T1_0,beta_T1(1));
                    warning('Change T1 Calibration')
                end
                if  abs(T2_0-T2_0_def)>1e-3 || abs(B2_0_def- beta_T2(1)>1e-2)
                    temp = regexprep(temp,sprintf('(SN(.){0,11}= %5s(((.|\\n){0,100})%7.3f))',char(setupstr(raw_mat.cfgobj, 'T2','SN')),str2double(setupstr(raw_mat.cfgobj, 'T2','T_0'))),...
                        sprintf(['SN          = %5s\n'...
                        'beta_1      = %7.2f\n'...
                        'T_0         = %7.3f\n'],...
                        char(setupstr(raw_mat.cfgobj, 'T2','SN')),beta_T2(1),T2_0));
                    fprintf(logid,'%19s %8s T2 Calibration T_0: %7.3f  Beta1: %7.2f   \n',datestr(now,'yyyy/mm/dd HH:MM:SS'),Raw_list(i).name,T2_0,beta_T2(1));
                    warning('Change T2 Calibration')
                end
                
                cfgid = fopen('./setup_changT.cfg','w');
                fprintf(cfgid,temp);
                fclose(cfgid);
                patch_setupstr([VMP_RAWP_Path Raw_list(i).name],'setup_changT.cfg')
                delete('./setup_changT.cfg')
                %%% reprocess
                warning('Reprocess FP07 Calibration')
                raw_mat = odas_p2mat([VMP_RAWP_Path Raw_list(i).name]);
                if mean(raw_mat.T1_fast)<0
                    warning('T1 is broken and replaced with T2 (pause) ')
                    pause
                    raw_mat.T1 = raw_mat.T2;
                    raw_mat.T1_dT2 = raw_mat.T1_dT2;
                    raw_mat.T1_fast = raw_mat.T2_fast;
                    raw_mat.T1_slow = raw_mat.T2_slow;
                elseif mean(raw_mat.T2_fast)<0
                    warning('T2 is broken and replaced with T1 (pause)')
                    pause
                    raw_mat.T2 = raw_mat.T1;
                    raw_mat.T2_dT2 = raw_mat.T1_dT1;
                    raw_mat.T2_fast = raw_mat.T1_fast;
                    raw_mat.T2_slow = raw_mat.T1_slow;
                elseif mean(raw_mat.T1_fast)<0 && mean(raw_mat.T2_fast)<0
                    error('Two FP07 sensors are broken!!!')
                end
                
                save([VMP_RAWM_Path Raw_list(i).name(1:end-2) '.mat'],'-struct','raw_mat','-v7.3')
                logid = fopen('./vmp_process_log.txt','a');
                fprintf(logid,'%19s %8s Converting P-file.... \n',datestr(now,'yyyy/mm/dd HH:MM:SS'),Raw_list(i).name);
                
                cali = 1;
            else
                cali = 0;
            end
        end
        
        %% Jac salinity lag cali. & Density
        %match JAC_T w/ JAC_C&FP07
        lag_JAC = dsearchn(xcorr(raw_mat.T1_slow,raw_mat.JAC_T,20,'normalized'),1)-21;
        raw_mat.JAC_T = circshift(raw_mat.JAC_T,lag_JAC);
        
        raw_mat.JAC_SP = gsw_SP_from_C(raw_mat.JAC_C,raw_mat.JAC_T,raw_mat.P_slow);
        raw_mat.JAC_SA = gsw_SA_from_SP(raw_mat.JAC_SP,raw_mat.P_slow,lon_c,lat_c);
        raw_mat.JAC_theta = gsw_CT_from_t(raw_mat.JAC_SA,raw_mat.JAC_T,raw_mat.P_slow); %Conservation T
        raw_mat.JAC_sigma = gsw_sigma0(raw_mat.JAC_SA,raw_mat.JAC_theta);
        raw_mat.JAC_rho = gsw_rho(raw_mat.JAC_SA,raw_mat.JAC_theta,raw_mat.P_slow)-1000; %subtract 1000 kg/m^3
        
        %% print profile
        %     figure
        %     plot(raw_mat.t_slow,-raw_mat.P_slow,'linewidth',2)
        %     hold on
        %     for pro_i = 1:length(profile_idx)
        %         plot(raw_mat.t_slow(profile_idx(1,pro_i):profile_idx(2,pro_i)),-raw_mat.P_slow(profile_idx(1,pro_i):profile_idx(2,pro_i)),'r','linewidth',2)
        %     end
        %     hold off
        %     xlabel('time')
        %     ylabel('pressure')
        %     pause
        %     close all
        
        %% calculating dissipation rate & create file
        % For more than 2 SH/T sensor OR C-sensor OR diff vehicle please refer
        % quick_look to alter the code
        
        [slow_matrix,fast_matrix] = deal(nan(max(diff(profile_idx)),size(profile_idx,2)));
        
        for j = 1:size(profile_idx,2)
            pr_str = profile_idx(1,j); %profile start point
            
            pr_end = profile_idx(2,j); %profile end point
            pr_idx = (pr_str:pr_end);
            
            [~,pr_str_fast] = min(abs(raw_mat.t_fast-raw_mat.t_slow(pr_str)));
            [~,pr_end_fast] = min(abs(raw_mat.t_fast-raw_mat.t_slow(pr_end)));
            pr_idx_fast = (pr_str_fast:pr_end_fast);
            
            if ~exist([VMP_PROC_Path Prefix '_' datestr(raw_mat.filetime+raw_mat.t_slow(pr_str)/86400,'yyyymmddHHMMSS') '.mat'],'file')
                %% despike
                % piezo-accelerometers
                [Ax_des, spikes_Ax, pass_count_Ax, raction_Ax] = despike(raw_mat.Ax(pr_idx_fast), despike_A.thres, despike_A.smooth , raw_mat.fs_fast, round(despike_A.N_FS*raw_mat.fs_fast));
                [Ay_des, spikes_Ay, pass_count_Ay, raction_Ay] = despike(raw_mat.Ay(pr_idx_fast), despike_A.thres, despike_A.smooth , raw_mat.fs_fast, round(despike_A.N_FS*raw_mat.fs_fast));
                
                if raction_Ax*100 >2|| raction_Ay*100>2
                    warning('Spike ratio: Ax: %5.1f %% Ay:%5.1f %% (Dat:%03i;cast:%03i)\n',raction_Ax*100,raction_Ay*100,Raw_list(i).name,j)
                    fprintf(logid,'%19s %8s Cast:%03i Piezo Spiking Ax: %5.1f %% Ay:%5.1f %% \n',datestr(now,'yyyy/mm/dd HH:MM:SS'),Raw_list(i).name,j,raction_Ax*100,raction_Ay*100);
                end
                
                AA = [Ax_des,Ay_des]; % accelerometer metrix for get_diss_odas
                
                % shear probe
                [SH1_des, spikes_sh1, pass_count_sh1, raction_sh1] = despike(raw_mat.sh1(pr_idx_fast), despike_sh.thres, despike_sh.smooth , raw_mat.fs_fast, round(despike_sh.N_FS*raw_mat.fs_fast));
                [SH2_des, spikes_sh2, pass_count_sh2, raction_sh2] = despike(raw_mat.sh2(pr_idx_fast), despike_sh.thres, despike_sh.smooth , raw_mat.fs_fast, round(despike_sh.N_FS*raw_mat.fs_fast));
                
                if raction_sh2*100 >3 || raction_sh1*100>3
                    warning('Spike ratio: SH1: %5.1f %% SH2:%5.1f %% (Dat:%03i;cast:%03i)\n',raction_sh1*100,raction_sh2*100,Raw_list(i).name,j)
                    fprintf(logid,'%19s %8s Cast:%03i Piezo Spiking Ax: %5.1f %% Ay:%5.1f %% \n',datestr(now,'yyyy/mm/dd HH:MM:SS'),Raw_list(i).name,j,raction_sh1*100,raction_sh2*100);
                end
                %% high_pass signal for dissipation
                HP_cut = 0.5 * 1/diss_info.fft_length_sec; %% follow Matlab Manual.
                [bh,ah] = butter(1, HP_cut/(raw_mat.fs_fast/2), 'high');
                
                %SH1
                SH1_HP = filter(bh, ah, SH1_des);
                SH1_HP = flipud(SH1_HP);
                SH1_HP = filter(bh, ah, SH1_HP);
                SH1_HP = flipud(SH1_HP);
                %SH2
                SH2_HP = filter(bh, ah, SH2_des);
                SH2_HP = flipud(SH2_HP);
                SH2_HP = filter(bh, ah, SH2_HP);
                SH2_HP = flipud(SH2_HP);
                
                %% dissipation
                diss_info.fft_length    = round(diss_info.fft_length_sec.*raw_mat.fs_fast);
                diss_info.diss_length  = diss_info.diss_length_fac.*diss_info.fft_length;
                diss_info.overlap = ceil(diss_info.diss_length/2);
                diss_info.fs_fast = raw_mat.fs_fast;
                diss_info.fs_slow = raw_mat.fs_slow;
                diss_info.speed = raw_mat.speed_fast(pr_idx_fast);
                diss_info.T = (raw_mat.T1_fast(pr_idx_fast)*2+0*raw_mat.T2_fast(pr_idx_fast))/2;
                diss_info.t = raw_mat.t_fast(pr_idx_fast);
                diss_info.P = raw_mat.P_fast(pr_idx_fast);
                
                diss = get_diss_odas([SH1_HP,SH2_HP],AA,diss_info);
                
                diss.P_diss = diss.P;
                diss = rmfield(diss,'P');
                
                %% diss check
                temp = length(find(abs(log10(diss.e(1,:)./diss.e(2,:)))>log10(5)))/length(diss.e(1,:));
                if temp>0.15
                    warning('Bad dissipation data %5.1f %% (Dat:%03i;cast:%03i)',temp*100,i,j)
                    fprintf(logid,'%19s %8s Cast:%03i disp ratio>5: %5.1f %% \n',datestr(now,'yyyy/mm/dd HH:MM:SS'),Raw_list(i).name,j,temp*100);
                end
                
                %% bio data
                if isfield(raw_mat,'Chlorophyll')
                    bio.Chlorophyll = raw_mat.Chlorophyll(pr_idx_fast);
                    bio.Turbidity = raw_mat.Turbidity(pr_idx_fast);
                else
                    bio.Chlorophyll = nan(size(raw_mat.P_fast(pr_idx_fast)));
                    bio.Turbidity = nan(size(raw_mat.P_fast(pr_idx_fast)));
                end
                
                %% bottom stress module
                if bot_mod
                    VMP_BBL_module
                end
                
                %% attach GPS
                vmp.dn = raw_mat.filetime + raw_mat.t_slow(pr_str)/86400;
                
                [~,time_idx] = min(abs(ship.dn-vmp.dn));
                vmp.lat = ship.lat(time_idx,1);
                vmp.lon = ship.lon(time_idx,1);
                vmp.dist_vmp = ship.dist_ship(time_idx,1);
                
                %% prepare mat-file
                % name save into vmp mat
                variable_names=[{'Ax'};{'Ay'};{'JAC_theta'};{'JAC_C'};{'JAC_SP'};{'JAC_SA'};{'JAC_sigma'};{'JAC_rho'};{'JAC_T'};...
                    {'T1_fast'};{'T2_fast'};{'T1_slow'};{'T2_slow'};{'gradT1'};{'gradT2'},;...
                    {'W_slow'};{'P_slow'};{'P_fast'};{'Incl_X'};{'Incl_Y'}];
                
                for k=1:length(variable_names)
                    eval(['N=length(raw_mat.' char(variable_names(k)) ');'])
                    if N==length(raw_mat.P_slow)
                        eval(['vmp.' char(variable_names(k)) '=raw_mat.' char(variable_names(k)) '(pr_idx);'])
                    elseif N==length(raw_mat.P_fast)
                        eval(['vmp.' char(variable_names(k)) '=raw_mat.' char(variable_names(k)) '(pr_idx_fast);'])
                    end
                end
                vmp.Ax = Ax_des;
                vmp.Ay = Ay_des;
                vmp.SH1 = SH1_des;
                vmp.SH2 = SH2_des;
                vmp.sn = str2double(setupstr(raw_mat.cfgobj,'instrument_info','sn'));
                vmp.data = raw_mat.input_parameters.fname;
                vmp.data_num = str2double(raw_mat.input_parameters.fname(end-4:end-2)); %data # on P-file; ex Dat_#.p
                vmp.profile = j; % profile number during each deployment
                
                %% save mat-file
                
                save([VMP_PROC_Path Prefix '_' datestr(raw_mat.filetime+raw_mat.t_slow(pr_str)/86400,'yyyymmddHHMMSS') '.mat'],'vmp','diss_info','diss','bio')
                
                if bot_mod
                    save([VMP_PROC_Path Prefix '_' datestr(raw_mat.filetime+raw_mat.t_slow(pr_str)/86400,'yyyymmddHHMMSS') '.mat'],'BBL_info','BBL','-append')
                end
            else
                load([VMP_PROC_Path Prefix '_' datestr(raw_mat.filetime+raw_mat.t_slow(pr_str)/86400,'yyyymmddHHMMSS') '.mat'])
            end
            %% write combined file
            %% VMP JAC&T1/T2
            vmp_names=[{'JAC_theta'};{'JAC_T'};{'JAC_C'};{'JAC_SP'};{'JAC_SA'};{'JAC_sigma'};{'JAC_rho'};...
                {'T1_fast'};{'T2_fast'};{'T1_slow'};{'T2_slow'};{'gradT1'};{'gradT2'},;...
                {'W_slow'};{'P_slow'};{'P_fast'}];
            
            N_slow = length(pr_idx);
            N_fast = length(pr_idx_fast);
            
            for k = 1:length(vmp_names)
                N = length(vmp.(vmp_names{k}));
                if j == 1
                    if N==N_slow
                        vmp_profile.(vmp_names{k}) = slow_matrix;
                        vmp_profile.(vmp_names{k})(1:N_slow,j) = vmp.(vmp_names{k});
                    elseif N==N_fast
                        vmp_profile.(vmp_names{k}) = fast_matrix;
                        if N_fast>size(vmp_profile.(vmp_names{k}),1)
                            vmp_profile.(vmp_names{k})(size(vmp_profile.(vmp_names{k}),1)+1:N_fast,:) = nan;
                        end
                        vmp_profile.(vmp_names{k})(1:N_fast,j) = vmp.(vmp_names{k});
                    end
                else
                    if N>size(vmp_profile.(vmp_names{k}),1) && N==N_fast
                        vmp_profile.(vmp_names{k})(size(vmp_profile.(vmp_names{k}),1)+1:N_fast,:) = nan;
                    end
                    vmp_profile.(vmp_names{k})(1:N,j) = vmp.(vmp_names{k});
                    vmp_profile.(vmp_names{k})(vmp_profile.(vmp_names{k}) == 0) = nan;
                end
            end
            
            %% bio data
            bio_names = {'Chlorophyll','Turbidity'};
            
            for k = 1:length(bio_names)
                if j == 1
                    vmp_profile.(bio_names{k}) = fast_matrix;
                    vmp_profile.(bio_names{k})(1:N_fast,j) = bio.(bio_names{k});
                else
                    if N_fast>size(vmp_profile.(vmp_names{k}),1)
                        vmp_profile.(vmp_names{k})(size(vmp_profile.(vmp_names{k}),1)+1:N_fast,:) = nan;
                    end
                    vmp_profile.(bio_names{k})(1:N_fast,j) = bio.(bio_names{k});
                    vmp_profile.(bio_names{k})(vmp_profile.(bio_names{k}) == 0) = nan;
                end
            end
            
            %% dissipation rate & figure of matrix
            diss_names = [{'e'},{'FM'}];
            N_diss = length(diss.P_diss);
            for k = 1:length(diss_names)
                if j == 1
                    vmp_profile.(diss_names{k}) = nan(2,length(diss.P_diss),size(profile_idx,2));
                    vmp_profile.(diss_names{k})(:,:,j) = diss.(diss_names{k});
                else
                    if  N_diss>size(vmp_profile.(diss_names{k}),2)
                        vmp_profile.(diss_names{k})(:,size(vmp_profile.(diss_names{k}),2)+1:N_diss,:) = nan;
                    end
                    vmp_profile.(diss_names{k})(:,1:N_diss,j) = diss.(diss_names{k});
                    vmp_profile.(diss_names{k})(vmp_profile.(diss_names{k}) == 0) = nan;
                end
            end
            
            diss_names = {'P_diss'};
            N_diss = length(diss.P_diss);
            for k = 1:length(diss_names)
                if j == 1
                    vmp_profile.(diss_names{k}) = nan(length(diss.P_diss),size(profile_idx,2));
                    vmp_profile.(diss_names{k})(:,j) = diss.(diss_names{k});
                else
                    if  N_diss>size(vmp_profile.(diss_names{k}),1)
                        vmp_profile.(diss_names{k})(size(vmp_profile.(diss_names{k}),2)+1:N_diss,:) = nan;
                    end
                    vmp_profile.(diss_names{k})(1:N_diss,j) = diss.(diss_names{k});
                    vmp_profile.(diss_names{k})(vmp_profile.(diss_names{k}) == 0) = nan;
                end
            end
            
            %% time/lat/lon/dist/data#
            vmp_profile.sn(j) = vmp.sn;
            vmp_profile.data_num(j) = vmp.data_num; %data # on P-file; ex Dat_#.p
            vmp_profile.profile(j) = vmp.profile; % profile number during each deployment
            
            vmp_profile.dn(j) = vmp.dn;
            vmp_profile.lat(j) = vmp.lat;
            vmp_profile.lon(j) = vmp.lon;
            vmp_profile.dist_vmp(j) = vmp.dist_vmp;
            
            if isstruct(BBL_info)
                vmp_profile.u_star(j) = BBL.u_star;
                vmp_profile.u_star_cint(j,:) = BBL.u_star_cint;
                vmp_profile.tau(j) = BBL.tau;
                vmp_profile.ToB(j) = BBL.ToB;
            end
            
        end
        save([VMP_PROC_P_Combine_Path  Prefix '_raw_profile_combine' datestr(vmp_profile.dn(end),'yyyymmddHHMMSS') '.mat'],'vmp_profile')
        
        %% bin
        vmp_names=[{'JAC_theta'};{'JAC_T'};{'JAC_C'};{'JAC_SP'};{'JAC_SA'};{'JAC_sigma'};{'JAC_rho'};{'e'};{'Chlorophyll'};{'Turbidity'}];
        vmp_names_combo=[{'theta'};{'T'};{'C'};{'SP'};{'SA'};{'sigma'};{'rho'};{'epsi'};{'Fl'};{'Turbi'}]; %%% Unified Name for Project
        
        for k = 1:length(vmp_names)
            temp = vmp_profile.(vmp_names{k});
            
            if size(vmp_profile.(vmp_names{k}),1) == size(vmp_profile.P_slow,1)
                temp3 = nan(level,size(vmp_profile.P_slow,2));
                for z = 1:level
                    mask = nan(size(vmp_profile.P_slow));
                    mask(vmp_profile.P_slow >= (z-1)*dz & vmp_profile.P_slow < z*dz) = 1;
                    temp3(z,:) = nanmean(mask.*temp);
                end
                temp3(temp3==0) = nan;
                
            elseif size(temp,1) == size(vmp_profile.P_fast,1)
                temp3 = nan(level,size(vmp_profile.P_fast,2));
                for z = 1:level
                    mask = nan(size(vmp_profile.P_fast));
                    mask(vmp_profile.P_fast >= (z-1)*dz & vmp_profile.P_fast < z*dz) = 1;
                    temp3(z,:) = nanmean(mask.*temp);
                end
                temp3(temp3==0) = nan;
                
            elseif size(temp,2) == size(vmp_profile.P_diss,1)
                mask = squeeze(temp(1,:,:)./temp(2,:,:)>5 | temp(1,:,:)./temp(2,:,:)<1/5);
                temp2 = squeeze(mean(temp,1));
                temp4 = min(temp(1,:,:),temp(2,:,:));
                temp2(mask) = temp4(mask);
                
                temp3 = nan(level,size(vmp_profile.P_diss,2));
                for z = 1:level
                    mask = nan(size(vmp_profile.P_diss));
                    mask(vmp_profile.P_diss >= (z-1)*dz & vmp_profile.P_diss < z*dz) = 1;
                    temp3(z,:) = nanmean(mask.*temp2);
                end
                temp3(temp3==0) = nan;
            end
            vmp_combo_temp.(vmp_names_combo{k})= temp3;
        end
        
        vmp_combo_temp.dn = vmp_profile.dn;
        vmp_combo_temp.lat = vmp_profile.lat;
        vmp_combo_temp.lon = vmp_profile.lon;
        vmp_combo_temp.dist_vmp = vmp_profile.dist_vmp;
        vmp_combo_temp.data_num = vmp_profile.data_num;
        vmp_combo_temp.profile = vmp_profile.profile;
        vmp_combo_temp.depth = (-0.5:-1:-level+0.5);
        
        if isfield(vmp_profile,'tau')
            vmp_combo_temp.u_star = vmp_profile.u_star;
            vmp_combo_temp.u_star_cint = vmp_profile.u_star_cint;
            vmp_combo_temp.tau = vmp_profile.tau;
            vmp_combo_temp.ToB = vmp_profile.ToB;
        end
        
        
        
        if exist([VMP_PROC_final_Path Prefix 'VMP_Precess.mat'],'file')
            vmp_combo = load([VMP_PROC_final_Path Prefix 'VMP_Precess.mat']);
            
            N_str = length(vmp_combo.dn)+1;
            
            vmp_names_combo=[{'theta'};{'T'};{'C'};{'SP'};{'SA'};{'sigma'};{'rho'};{'epsi'};{'Fl'};{'Turbi'}]; %%% Unified Name for Project
            
            for k = 1:length(vmp_names_combo)
                vmp_combo.(vmp_names_combo{k})(:,N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.(vmp_names_combo{k});
            end
            
            vmp_combo.dn(N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.dn;
            vmp_combo.lat(N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.lat;
            vmp_combo.lon(N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.lon;
            vmp_combo.dist_vmp(N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.dist_vmp;
            vmp_combo.data_num(N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.data_num;
            vmp_combo.profile(N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.profile;
            vmp_combo.depth = (-dz/2:-dz:-dz*level+dz/2);
            
            if isfield(vmp_profile,'tau')
                vmp_combo.u_star(N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.u_star;
                vmp_combo.u_star_cint(N_str:N_str-1+length(vmp_combo_temp.dn),:) = vmp_combo_temp.u_star_cint;
                vmp_combo.tau(N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.tau;
                vmp_combo.ToB(N_str:N_str-1+length(vmp_combo_temp.dn)) = vmp_combo_temp.ToB;
            end
        else
            vmp_combo = vmp_combo_temp;
        end
        save([VMP_PROC_final_Path Prefix '_VMP_Precessed.mat'],'-struct','vmp_combo','-v7.3')
        movefile([VMP_RAWP_Path Raw_list(i).name],[VMP_RAWP_Path '/done/' Raw_list(i).name])
        clear vmp_profile vmp_combo_temp vmp_combo
    end
end
fclose(logid);

%% sort and remove duplicated profile

vmp_combo = load([VMP_PROC_final_Path Prefix '_VMP_Precessed.mat']);

[~,idx_temp] = sort(vmp_combo.dn);

duplicate_idx = find(diff(vmp_combo.dn(idx_temp))==0);

idx_temp(duplicate_idx) = [];

vmp_names_combo=[{'theta'};{'T'};{'C'};{'SP'};{'SA'};{'sigma'};{'rho'};{'epsi'};{'Fl'};{'Turbi'}]; %%% Unified Name for Project

for k = 1:length(vmp_names_combo)
    vmp_combo.(vmp_names_combo{k}) = vmp_combo.(vmp_names_combo{k})(:,idx_temp);
end

vmp_combo.dn = vmp_combo.dn(idx_temp);
vmp_combo.lat = vmp_combo.lat(idx_temp);
vmp_combo.lon = vmp_combo.lon(idx_temp);
vmp_combo.dist_vmp = vmp_combo.dist_vmp(idx_temp);
vmp_combo.data_num = vmp_combo.data_num(idx_temp);
vmp_combo.profile = vmp_combo.profile(idx_temp);

if isfield(vmp_combo,'tau')
    vmp_combo.u_star = vmp_combo.u_star(idx_temp);
    vmp_combo.u_star_cint= vmp_combo.u_star_cint(idx_temp,:);
    vmp_combo.tau = vmp_combo.tau(idx_temp);
    vmp_combo.ToB = vmp_combo.ToB(idx_temp);
end

save([VMP_PROC_final_Path Prefix '_VMP_Precessed.mat'],'-struct','vmp_combo','-v7.3')
