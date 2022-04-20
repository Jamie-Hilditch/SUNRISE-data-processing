function out = WWvel_downward_SR(WWmeta,variables,k)
% out{1} = [];  % E-W velocity
% out{2} = [];  % N-s velocity
% out{3} = [];  % time
% out{4} = [];  % depth
% bofu zheng

%% set variables
if nargin<2
    num      = 20;   % 20 files in a group 
    blockdis = 0.1;  % blocking distance 0.1 m
    cellsize = 0.25; % cell size 0.25 m
    saprate  = 16;   % sampling rate 16 Hz
    boxs  = 0.25; % box size 0.25 m
    z_max    = 100;  % profile depth 100 m
else
    num      = variables.NUM_combining_files;
    blockdis = variables.blockdis;
    cellsize = variables.cellsize;
    saprate  = variables.saprate;
    boxs  = variables.boxsize;
    z_max    = variables.z_max;
end

%%

out{1} = [];  % E-W velocity
out{2} = [];  % N-s velocity
out{3} = [];  % time
out{4} = [];  % depth

dd0 = dir([WWmeta.aqdpath '*.mat']);
inds = 1:num:length(dd0);  % starting index for each group
inde = inds(2:end)-1;
inde = [inde length(dd0)]; % ending index for each group


vele_up = [];
veln_up = [];
time_up = [];
for q = 1:length(inds)
    filename = ['Profiles_downcast_',WWmeta.name_aqd,'_',num2str(inds(q)),'_',num2str(inde(q)),'.mat']; % Name of the processed file1
    load([WWmeta.propath_rearrange filename]);
    % find long enough profile
%     AQDprofiles_up = AQDprofiles;
%     index=find(abs(cellfun(@(x) x.Burst_Pressure(end)-x.Burst_Pressure(1),AQDprofiles_up))>0.4*z_max);
    index=find(abs(cellfun(@(x) x.Burst_Pressure(end)-x.Burst_Pressure(1),AQDprofiles_down))>0.4*z_max);
    
    %% load data
    for k = 1:length(index)
%         time = AQDprofiles_up{index(k)}.Burst_Time;
%         dpth = AQDprofiles_up{index(k)}.Burst_Pressure;
%         corr1 = double(AQDprofiles_up{index(k)}.Burst_CorBeam1);
%         corr2 = double(AQDprofiles_up{index(k)}.Burst_CorBeam2);
%         corr3 = double(AQDprofiles_up{index(k)}.Burst_CorBeam3);
%         corr4 = double(AQDprofiles_up{index(k)}.Burst_CorBeam4);
%         vele{1} = AQDprofiles_up{index(k)}.Burst_VelEast;
%         vele{2} = AQDprofiles_up{index(k)}.Burst_VelNorth;
%         vele{3} = 1/2*(AQDprofiles_up{index(k)}.Burst_VelUp1 + ...
%             AQDprofiles_up{index(k)}.Burst_VelUp2);  % raw vel z
%         pitch = AQDprofiles_up{index(k)}.Burst_Pitch;
%         roll = AQDprofiles_up{index(k)}.Burst_Roll;
%         heading = AQDprofiles_up{index(k)}.Burst_Heading;
%         acce = AQDprofiles_up{index(k)}.Burst_Accelerometer;
        
        time = AQDprofiles_down{index(k)}.Burst_Time;
        dpth = AQDprofiles_down{index(k)}.Burst_Pressure;
        corr1 = double(AQDprofiles_down{index(k)}.Burst_CorBeam1);
        corr2 = double(AQDprofiles_down{index(k)}.Burst_CorBeam2);
        corr3 = double(AQDprofiles_down{index(k)}.Burst_CorBeam3);
        corr4 = double(AQDprofiles_down{index(k)}.Burst_CorBeam4);
        vele{1} = AQDprofiles_down{index(k)}.Burst_VelEast;
        vele{2} = AQDprofiles_down{index(k)}.Burst_VelNorth;
        vele{3} = 1/2*(AQDprofiles_down{index(k)}.Burst_VelUp1 + ...
            AQDprofiles_down{index(k)}.Burst_VelUp2);  % raw vel z
        pitch = AQDprofiles_down{index(k)}.Burst_Pitch;
        roll = AQDprofiles_down{index(k)}.Burst_Roll;
        heading = AQDprofiles_down{index(k)}.Burst_Heading;
        acce = AQDprofiles_down{index(k)}.Burst_Accelerometer;
     
                
        NCells = AQDprofiles_down{index(k)}.Burst_NCells(1); % number of cells
        Npings = length(time);   % number of pings
     
        % extract velocity that has low correlation
        m = corr1<40 | corr2<40 | corr3<40 | corr4<40;
        vele{1}(find(m == 1)) = nan;
        vele{2}(find(m == 1)) = nan;
        vele{3}(find(m == 1)) = nan;
       
        % calculate ww motion velocity
        velemc = WWcorr(time,dpth,acce,pitch,roll,heading,vele,saprate);
            
        % build dpth matrix
        dpth_temp = -dpth*ones(1,NCells)-ones(Npings,1)*cellsize*double(0:1:NCells-1) - blockdis;
        dpth_temp(find(m == 1)) = nan;
        dpth_temp_all = dpth_temp(:);
        
        % build time matrix
        time_temp = time*ones(1,NCells);
        time_temp(find(m==1)) = nan;
        time_temp_all = time_temp(:);
        
        %% new for SR
        %         [~, Id]=min(abs(dpth(:)-10));   % find index for depth at 7m
        %
        %         dpth_temp_2 = -dpth(1:Id)*ones(1,NCells)-ones(Id,1)*cellsize*double(0:1:NCells-1) - blockdis;
        %         dpth_temp_all_2 = dpth_temp_2(:);
        %
        %         ucon_all = velemc{1}(1:Id,:);  % to find constant velocity
        %         vcon_all = velemc{2}(1:Id,:);
        %
        %         index_box_2 = find(-20<=dpth_temp_all_2 & dpth_temp_all_2<=-18); % -20 and -17 previously
        %         ucon = nanmean(ucon_all(index_box_2));
        %         vcon = nanmean(vcon_all(index_box_2));
        %
        %         for i = 1:length(time)
        %             index_box_3 = find(-20<=dpth_temp(i,:) & dpth_temp(i,:)<=-18);  % should be small values
        %             uav = nanmean(velemc{1}(i,index_box_3));
        %             vav = nanmean(velemc{2}(i,index_box_3));
        %
        %             velemc{1}(i,:) = velemc{1}(i,:) - (uav-ucon);
        %             velemc{2}(i,:) = velemc{2}(i,:) - (vav-vcon);
        %         end
        
        %         velemc{1} = velemc{1}-ucon;
        %         velemc{2} = velemc{2}-vcon;
        
        
        %
        [~, Id]=min(abs(dpth(:)-10));   % find index for depth at 10m
        lid = length(dpth)-Id+1;
        dpth_temp_2 = -dpth(Id:end)*ones(1,NCells)-ones(lid,1)*cellsize*double(0:1:NCells-1) - blockdis;
        dpth_temp_all_2 = dpth_temp_2(:);
        
        ucon_all = velemc{1}(Id:end,:);  % to find constant velocity
        vcon_all = velemc{2}(Id:end,:);
        
        index_box_2 = find(-20<=dpth_temp_all_2 & dpth_temp_all_2<=-18);  % -20 - -17 previously
        ucon = nanmean(ucon_all(index_box_2));
        vcon = nanmean(vcon_all(index_box_2));
        
        for i = 1:length(time)
            index_box_3 = find(-20<=dpth_temp(i,:) & dpth_temp(i,:)<=-18);  % should be small values
            uav = nanmean(velemc{1}(i,index_box_3));
            vav = nanmean(velemc{2}(i,index_box_3));
            
            velemc{1}(i,:) = velemc{1}(i,:) - (uav-ucon);
            velemc{2}(i,:) = velemc{2}(i,:) - (vav-vcon);
        end

        %%
        % all velocity 
        velemcE_temp_all = velemc{1}(:);
        velemcN_temp_all = velemc{2}(:);
        
        dz = (0:-boxs:-z_max)';
        temp_pro = nan(length(dz),2);  % each profile, first column: E/W; second column: N/S
        temp_time   = nan(length(dz),1);  % box-averaged time
        for i = 1:length(dz)
            index_box = [];
            index_box = find(dz(i)-boxs/2<dpth_temp_all & dpth_temp_all<=dz(i)+boxs/2);
            if length(index_box)>0
                temp_pro(i,1)  = nanmean(velemcE_temp_all(index_box));
                temp_pro(i,2)  = nanmean(velemcN_temp_all(index_box));
                temp_time(i,1) = nanmean(time_temp_all(index_box));                       
            end
        end
        clear i
            
        vele_up = [vele_up temp_pro(:,1)];
        veln_up = [veln_up temp_pro(:,2)];
        time_up = [time_up temp_time];
        
        disp([filename,' total #= ',num2str(size(time_up,2))]);
    end
end

%%
out{1} = vele_up;
out{2} = veln_up;
out{3} = time_up;
out{4} = -dz*ones(1,length(time_up));


end


