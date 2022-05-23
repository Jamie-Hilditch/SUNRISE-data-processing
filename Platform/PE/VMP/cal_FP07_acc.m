function [T_0,beta,Lag] = cal_FP07_acc(file_mat,profile_idx,T_ref_string,T_string,order)
% This code is modified from RL's cal_FP07_in_situ.m. The code will
% compare more profiles than orignal code and enhence efficiency. It will
% provide more accurate calibration. The code is based on RL's 2021
% library, one should check if there is any update.
% The code is only for VMP!!!
% Fucent Hsu, CEOAS, OSU, 2/27/2022

T_with_pre_emphasis_string = [T_string '_d' T_string]; % Should be T1_dT1 or T2_dT2
T_without_pre_emphasis_string = T_string; % Should be T1 or T2


var_cell = {T_with_pre_emphasis_string,T_without_pre_emphasis_string,T_ref_string,'t_fast', 't_slow', 'P_slow', 'W_slow', 'fs_fast', 'fs_slow','setupfilestr'};
for i = 1:length(var_cell)
    eval([var_cell{i} '= file_mat.(var_cell{i});']);
end

% Assign temperature vectors
eval(['T_without_pre_emphasis = ' T_without_pre_emphasis_string ';'])

eval (['T_ref = ' T_ref_string ';']) % T_ref is the reference thermometer, usually SBT or JAC_T
eval(['T = ' T_without_pre_emphasis_string ';']) % T is the thermistor signal without pre-emphasis

% If we have a signal with pre-emphasis, then we will use it to form the signal T.
eval(['T = ' T_with_pre_emphasis_string ';']) % T is the thermistor signal with pre-emphasis
T = deconvolve(...
    T_with_pre_emphasis_string, ...
    T_without_pre_emphasis, ...
    T, ...
    fs_fast, setupfilestr);

% sampling rate ratio (used to downsample)
ratio = round(fs_fast / fs_slow);

% down size to match T_ref
T = reshape(T, ratio, []);
T = mean(T)';

m = false(size(T));
for i = 1:length(profile_idx)
    m(profile_idx(1,i):profile_idx(2,i)) = 1;
end
%----------------------------------
% - Warning if range is too small -
%----------------------------------
if max(T_ref(m))-min(T_ref(m))<=8 && order>1
    warning('read!!!')
    temp = input('Temperature range is less than 8 degrees. RL recommends using FIRST-ORDER calibration. \nDo you agree to change? [y,n]','s');
    if strcmp(temp,'y')
        order = 1;
    end
end

%-------------------------------------------------------------------------
% -- Low-Pass filter thermistor data and compute cross correlation -------
%-------------------------------------------------------------------------

% Low-pass filtering the thermistor data to make it more compatible with the JAC-T
if strcmp(T_ref_string,'JAC_T')
    W_mean = abs(mean(W_slow(m)));
    fc = 0.73 * sqrt(W_mean / 0.62); % in Hz
    [b,a] = butter(1, fc / (fs_slow/2));
    T = filter(b, a, T);
else
    fc = fs_slow/3; % It is a Sea-Bird Thermometer
    [b,a] = butter(1, fc / (fs_slow/2));
    T = filter(b, a, T);
end

% Compute cross-correlation
max_lag = round(10*fs_slow); % estimate of the max lag required to find the actual lag.
[bb, aa] = butter(2,4/(fs_slow/2)); % 4 Hz smoother to suppress high-frequency noise

for i = 1:length(profile_idx)
    [correlation, lags] = xcorr(...
        filter(bb,aa,detrend(diff(T(profile_idx(1,i):profile_idx(2,i))))),...
        filter(bb,aa,detrend(diff(T_ref(profile_idx(1,i):profile_idx(2,i))))),max_lag,'coeff');
    [max_corr, m_lag(i)] = max(abs(correlation));
end
m_lag = round(mean(m_lag));

m_lag = m_lag - max_lag - 1;
Lag    = m_lag / fs_slow; % in seconds and should be negative.

%-----------------------------------------------------
% -- Do regression to get thermistor coefficients ----
%           (using Steinhart-Hart equation)
%-----------------------------------------------------

% Copy only the profile data
T_ref_prof = T_ref(m);
T_prof     = T(m);
P_prof     = P_slow(m);

% First align the T and T_ref signals using m_lag.
if m_lag >0, m_lag = 0; end % m_lag is expected to be negative because
% reference sensor is physically 'behind' probes
P_prof     = P_prof(1:end+m_lag);
T_prof     = T_prof(1:end+m_lag);
T_ref_prof = T_ref_prof(1-m_lag:end); % shift reference temperature
T_ref_regress = T_ref_prof + 273.15; % in kelvin
T_ref_regress = 1 ./ T_ref_regress;

% Now gather information about the electronics for this thermistor.
therm_type =    (char(setupstr( file_mat.cfgobj, T_string, 'type')));
E_B = str2double(char(setupstr( file_mat.cfgobj, T_string, 'E_B')));
a   = str2double(char(setupstr( file_mat.cfgobj, T_string, 'a'  )));
b   = str2double(char(setupstr( file_mat.cfgobj, T_string, 'b'  )));
G   = str2double(char(setupstr( file_mat.cfgobj, T_string, 'G'  )));
adc_fs   = str2double(char(setupstr( file_mat.cfgobj, T_string, 'adc_fs'  )));
adc_bits = str2double(char(setupstr( file_mat.cfgobj, T_string, 'adc_bits'  )));
try zero = str2double(char(setupstr( file_mat.cfgobj, T_string, 'adc_zero'  )));catch, zero = 0; end

% Compute non-dimensional thermistor voltage
if strcmp(therm_type, 'therm')
    factor = (adc_fs / 2^adc_bits)*2 / (G*E_B);
    Z = factor*(T_prof - a)/b;
elseif strcmp(therm_type, 't_ms')
    Z = T_prof * (adc_fs/2^adc_bits) + zero;
    Z = ((Z - a)/b) *2 / (G*E_B);
end

% Compute resistance ratio for this thermistor.
RT_R0 = (1 - Z) ./ (1 + Z);
RT_R0 = log(RT_R0);

% Generate the coefficients for this thermistor.
beta = zeros(1,order);
p = polyfit(RT_R0, T_ref_regress, order);
pp = p; % save for later usage
p = 1 ./ p;
p = fliplr(p); % place in ascending order
T_0    = p(1);
for index = 2:order+1
    beta(index-1) = p(index);
end


end

