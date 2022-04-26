clear
close all

%% initial setup
clear
close all

%%% setting path and pre-loading files
addpath('../_Config')
Process_Mode = 'HydroCombo';
data_path %% all data path and library



%%%%%%% sort %%%%%%%

[~,idx] = sort([vmp_combo.dn ctd_combo.dn]);

%% theta
temp = [vmp_combo.theta ctd_combo.theta];
hydro_combo.theta = temp(:,idx);

%% SP
temp = [vmp_combo.SP ctd_combo.SP];
hydro_combo.SP = temp(:,idx);

%% SA
temp = [vmp_combo.SA ctd_combo.SA];
hydro_combo.SA = temp(:,idx);

%% sigma
temp = [vmp_combo.sigma ctd_combo.sigma];
hydro_combo.sigma = temp(:,idx);

%%  epsi
temp = [vmp_combo.epsi nan(size(ctd_combo.SA))];
hydro_combo.epsi = temp(:,idx);

%%  Chlorophyll
temp = [vmp_combo.Fl nan(size(ctd_combo.SA))];
hydro_combo.Fl = temp(:,idx);

%%  Turbidity
temp = [vmp_combo.Turbi nan(size(ctd_combo.SA))];
hydro_combo.Turbi = temp(:,idx);

%%  O2A
temp = [nan(size(vmp_combo.theta)) ctd_combo.DO2A];
hydro_combo.DO2A = temp(:,idx);

%%  O2R
temp = [nan(size(vmp_combo.theta)) ctd_combo.DO2R];
hydro_combo.DO2R = temp(:,idx);

if isfield(vmp_combo,'tau')
    %%  ustar
    temp = [vmp_combo.u_star nan(size(ctd_combo.lat))];
    hydro_combo.u_star = temp(:,idx);

    %%  u_star_cint
    temp = [vmp_combo.u_star_cint nan(2,length(ctd_combo.lat))];
    hydro_combo.u_star_cint = temp(:,idx);

    %%  tau
    temp = [vmp_combo.tau nan(size(ctd_combo.lat))];
    hydro_combo.tau = temp(:,idx);

    %%  tau
    temp = [vmp_combo.ToB nan(size(ctd_combo.lat))];
    hydro_combo.ToB = temp(:,idx);    
end

%% lat 
temp = [vmp_combo.lat ctd_combo.lat];
hydro_combo.lat = temp(:,idx);

%% lon 
temp = [vmp_combo.lon ctd_combo.lon];
hydro_combo.lon = temp(:,idx);

%% dn 
temp = [vmp_combo.dn ctd_combo.dn];
hydro_combo.dn = temp(:,idx);

%% dist 
temp = [vmp_combo.dist_vmp ctd_combo.dist_ctd];
hydro_combo.dist = temp(:,idx);

%% depth 
if ctd_combo.depth~=vmp_combo.depth
    error('CTD & VMP have different bin')
else
    hydro_combo.depth = vmp_combo.depth;
end


%% 
save([Hydro_DATA_Path Prefix '_HydroCombo_Processed.mat'],'-struct','hydro_combo','-v7.3')
