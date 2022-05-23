function out = WWcorr(time,dpth,acce,pitch,roll,heading,vele,saprate)
% function to calculate WW motion corrected velocity for down-ward looking
% ADCP on the WW
% saprate - sampling rate
% 
% out{1} - after motion correction (E-W)
% out{2} - after motion correction (N-S)
% out{3} - after motion correction (U-D)
% out{4} - ww motion (E-W)
% out{5} - ww motion (N-S)
% out{6} - ww motion (U-D)

% Bofu Zheng
% modified on Nov 2 2020
%% parameter
NCells = size(vele{1},2);  % number of cells

%% band pass acceleration
% static acceleration
ax = sin(pitch*pi/180);  % x component tilt
ay = sin(roll*pi/180).*cos(pitch*pi/180);   % y component tilt
az = cos(roll*pi/180).*cos(pitch*pi/180);   % z component tilt

% dynamic accelerometer
accex = acce(:,1) - ax;   % get rid of tilt
accey = acce(:,2) - ay;
accez = acce(:,3) - az;

% transfer acceleration from XYZ into ENU system
aENU = XYZ2ENU(accex, accey, accez, pitch, roll, heading);

% filter acceleration
aENUf = [];  % creat a fake time series to get rid of the edge effect
aENUf(:,1) = [flip(aENU(:,1)); aENU(:,1); flip(aENU(:,1))];
aENUf(:,2) = [flip(aENU(:,2)); aENU(:,2); flip(aENU(:,2))];
aENUf(:,3) = [flip(aENU(:,3)); aENU(:,3); flip(aENU(:,3))];

[b,a] = butter(1,[0.1 1.2]/(saprate/2),'bandpass');  % set up filter, 0.1-1.2 Hz
acu = filtfilt(b,a,double(aENUf(:,1)));
acu = detrend(acu,'constant');  % get rid of undesired trend
acu = acu(length(acu)*1/3+1:length(acu)*2/3);  % chunk back to original length
acu = acu - nanmean(acu);

acv = filtfilt(b,a,double(aENUf(:,2)));
acv = detrend(acv,'constant');
acv = acv(length(acv)*1/3+1:length(acv)*2/3);
acv = acv - nanmean(acv);

acw = filtfilt(b,a,double(aENUf(:,3)));
acw = detrend(acw,'constant');
acw = acw(length(acw)*1/3+1:length(acw)*2/3);
acw = acw - nanmean(acw);

%% calculate WW motion - translational velocity
WWu = zeros(size(time));
WWv = zeros(size(time));
WWw = zeros(size(time));

for i = 1:length(WWu)-1
    WWu(i+1) = WWu(i) + acu(i)*9.81*86400*(time(i+1)-time(i));
    WWv(i+1) = WWv(i) + acv(i)*9.81*86400*(time(i+1)-time(i));
    WWw(i+1) = WWw(i) + acw(i)*9.81*86400*(time(i+1)-time(i));
end

% WWu = detrend(WWu);
% WWv = detrend(WWv);
% WWw = detrend(WWw);

%% motion corrected velocity 
un_mc = vele{1}+WWu*ones(1,NCells);     % E-W
vn_mc = vele{2}+WWv*ones(1,NCells);     % N-S

dpdt = gradient(dpth,time*86400);       % WW vertical velocity
wn_mc = vele{3}-dpdt*ones(1,NCells);    % U-D

%%
out{1} = un_mc;  % motion corrected
out{2} = vn_mc;
out{3} = wn_mc;
out{4} = WWu;    % WW motion
out{5} = WWv;
out{6} = WWw; %-dpdt;

end
        