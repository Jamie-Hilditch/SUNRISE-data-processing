%%% Module for VMP
%%% Bottom Boundary Layer: stress (abs) & depth
%%% Perlin et al (2005)

BBL_info.vers = 'v1.0'; % 1/5/2022

%% Parameter
BBL_info.fft_length_sec    = 0.2; % unit second, falling rate~1.2m/s
BBL_info.diss_length_fac   = 2; % no unit, >2*fft_length Rec:>3
BBL_info.P = raw_mat.P_fast(pr_idx_fast);
BBL_info.dr_bbl = 0.01; % boundary layer sigma difference 
BBL_info.fc_sigma_bbl = 0.5/0.5; %0.5/period_cut(sec)

%% Detect Bottom
% when VMP hits bottom,accelerometer detects a strong spike which will be removed.
% moving mean/diff filter are applied to find the blank zone.
idx_AA1 = find(movmean(movmean(diff(AA(:,1)),20),5)==0);
idx_AA2 = find(movmean(movmean(diff(AA(:,2)),20),5)==0);

idx_AA1((BBL_info.P(idx_AA1)-max(BBL_info.P))<-1) = [];
idx_AA2((BBL_info.P(idx_AA2)-max(BBL_info.P))<-1) = [];

pr_end_fast_BBL =  min(min(idx_AA1),min(idx_AA2))-12; %%12 is the lag due to filters. 
if isempty(pr_end_fast_BBL)
    pr_end_fast_BBL = length(BBL_info.P);
end

clear idx_AA2 idx_AA1

%% remove surface
% force RL code to process whole near bottom signal
[~,pr_str_fast_BBL] = min(abs(diss_info.t - diss_info.t(1) -  mod(diss_info.t(pr_end_fast_BBL) - diss_info.t(1),BBL_info.fft_length_sec)));

pr_idx_fast_BBL = (pr_str_fast_BBL:pr_end_fast_BBL);

%% BBL Profile
%high_pass signal for dissipation
HP_sigma_BBL = 0.5 * 1/BBL_info.fft_length_sec; %% follow Matlab Manual.
[bh_BBL,ah_BBL] = butter(1, HP_sigma_BBL/(raw_mat.fs_fast/2), 'high');

%SH1
SH1_HP_BBL = flipud(filter(bh_BBL, ah_BBL,flipud(filter(bh_BBL, ah_BBL, SH1_des(pr_idx_fast_BBL)))));

%SH2
SH2_HP_BBL = flipud(filter(bh_BBL, ah_BBL,flipud(filter(bh_BBL, ah_BBL, SH2_des(pr_idx_fast_BBL)))));


AA_BBL = [Ax_des(pr_idx_fast_BBL),Ay_des(pr_idx_fast_BBL)];

BBL_info.fft_length   = round(BBL_info.fft_length_sec.*raw_mat.fs_fast);
BBL_info.diss_length  = BBL_info.diss_length_fac.*BBL_info.fft_length;
BBL_info.overlap      = 0;
BBL_info.fs_fast = raw_mat.fs_fast;
BBL_info.fs_slow = raw_mat.fs_slow;
BBL_info.speed = diss_info.speed(pr_idx_fast_BBL);
BBL_info.T = diss_info.T(pr_idx_fast_BBL);
BBL_info.t = diss_info.t(pr_idx_fast_BBL);
BBL_info.P = diss_info.P(pr_idx_fast_BBL);

BBL.diss = get_diss_odas([SH1_HP_BBL,SH2_HP_BBL],AA_BBL,BBL_info);

%% changing vertical corodinate: height from ground

temp = pr_idx(end);
BP_BBL = raw_mat.P_slow(temp); % bottom pressure
BBL.HfG = -raw_mat.P_slow(pr_idx)+BP_BBL; % slow 
BBL.HfG_e = -BBL.diss.P+BP_BBL;
%% smoothing sigma profile
[bh_BBL, ah_BBL] = butter(1, BBL_info.fc_sigma_bbl/(raw_mat.fs_slow/2), 'high');

BBL.sigma_lp_bbl = raw_mat.JAC_sigma(pr_idx) - flipud(filter(bh_BBL,ah_BBL,flipud(filter(bh_BBL,ah_BBL,raw_mat.JAC_sigma(pr_idx))))); % low-pass

sigma_ToB = mean(BBL.sigma_lp_bbl(BBL.HfG>=0.2 & BBL.HfG<=0.3))-BBL_info.dr_bbl; %20-30 cm mean minus boundary def; top of boundary layer

BBL.ToB = BBL.HfG(find((BBL.sigma_lp_bbl(1:end-1)-sigma_ToB).*(BBL.sigma_lp_bbl(2:end)-sigma_ToB)<0,1,'last')); % boundary layer depth

%% fitting epsilon profile for tau
% follow Perlin 2005 stratification is not significant 0.6 Z. 
temp = find(BBL.HfG_e/BBL.ToB<=0.6);

if length(temp)>1
    u_star_fit = fit(0.4*BBL.HfG_e(temp)/BBL.ToB,BBL.diss.e(temp),'a/x','DiffMinChange',1e-14,'TolFun',1e-15,'TolX',1e-15,'StartPoint',1e-10,'Lower',0); % 0.4 ~ von Ka ́rma ́n’s constant

    BBL.u_star = (coeffvalues(u_star_fit))^(1/3); % fit value
    BBL.u_star_cint = confint(u_star_fit).^(1/3); % confidence interval

elseif length(temp)==1
    BBL.u_star = (BBL.diss.e(temp)*0.4*BBL.HfG_e(temp)/BBL.ToB)^(1/3);
    BBL.u_star_cint = [nan,nan];
else
    BBL.u_star = nan;
    BBL.u_star_cint = [nan,nan];
end

BBL.tau = (1000+sigma_ToB+BBL_info.dr_bbl)*BBL.u_star^2;

