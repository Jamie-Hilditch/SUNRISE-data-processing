%collect UHDAS ADCP

%% initial setup
clear
close all

%%%
ADCP_Project_name = 'WS21163_Hetland_merge';
ADCP_name = {'wh1200','wh600','wh600_no_beam2'};

%%% setting path and pre-loading files
addpath('../_Config')
Process_Mode = 'ADCP_UHDAS';
data_path %% all data path and library

%%
for i = 1:length(ADCP_name)
    
    ADCP_temp = load_getmat(fullfile( [ADCP_PROC_Path ADCP_name{i} '/contour/'], 'allbins_'));
    
    time_temp = num2str(ADCP_temp.time');
    ADCP_temp.dn = datenum(time_temp,'yyyy  mm  dd  HH  MM  SS');
    
    %true_headding = repmat(interp1(ship.dn,ship.true_headding,ADCP_temp.dn),[1 size(ADCP_temp.depth,1)])';
    
    %ADCP_temp.uA =  ADCP_temp.u.*cosd(true_headding) + ADCP_temp.v.*sind(true_headding);
    %ADCP_temp.uC =  ADCP_temp.u.*sind(true_headding) - ADCP_temp.v.*cosd(true_headding);    
    
    save([ADCP_PROC_final_Path Prefix '_' ADCP_name{i} '_Processed.mat'],'-struct','ADCP_temp','-v7.3')
    
    
end