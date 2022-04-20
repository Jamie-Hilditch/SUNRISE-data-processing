function [up,down,dataup,datadown] = mod_telemetry_getcastctd(CTD,min_depth,crit_depth)

% extract upcast downcast.
%  we filt the pressure with low pass at 1min. then look for the period
%  where pressure is higher ther min_depth (user defined).The downcast is
%  the part before the maximun pressure whithin this period the upscast is
%  the remaining. This is done iteratively until the end of the time serie. 
%  At the end, we only keep the profile where
%  max_depth-min_depth>=crit_depth. By default crit_depth is 10 m

% output: 
%   up:       upcasts indexes from the ctd
%   down:     downcasts indexes from the ctd
%   dataup:   data (P,S,T,C,sig,ctdtime) for the upcasts
%   datadown: data (P,S,T,C,sig,ctdtime) for the downcasts

%  Created by Arnaud Le Boyer on 7/28/18.

if nargin<3
    crit_depth=10;
end
% load(fullfile(CTD.CTDpath,CTD.name_ctd),CTD.name_ctd);
% eval(sprintf('data=%s;',CTD.name_ctd));
data = CTD;
% rename variable to make it easy, only for epsiWW
if isfield(data,'info')
    info=data.info;
    data=rmfield(data,'info');
end

if isfield(data,'sn')
    sn=data.sn;
    data=rmfield(data,'sn');
end

if isfield(data,'download')
    download=data.download;
    data=rmfield(data,'download');
end


% if statement to handle EPSI and WW processing
if ~isfield(data,'ctdtime')
    tdata=data.time;
else
    tdata=data.ctdtime;
end

% get time resolution
dt=median(diff(tdata)); % sampling period

% check data length
T=tdata(end)-tdata(1);
disp('check if time series is shorter than 3 hours')
if T<3/24  
    warning('time serie is less than 3 hours, very short for data processing, watch out the results')
end

% inverse pressure and remove outliers above the std deviation of a sliding window of 20 points
% pdata=filloutliers(data.P,'center','movmedian',20);
pdata = double(data.P);
if min(pdata)>8
    warning('there is an offset on the RBR pressure. We substract 8 m to get the profile.');
end
% buid a filter
disp('smooth the pressure to define up and down cast')
Nb  = 3; % filter order
fnb = 1/(2*dt); % Nyquist frequency
fc  = 1/600/dt; % 50 dt (give "large scale patern") 
[b,a]= butter(Nb,fc/fnb,'low');
% filt pdata
filt_pdata=filtfilt(b,a,pdata);
dfilt_pdata=diff(filt_pdata);
pddfilt_pdata = dfilt_pdata(1:end-1).*dfilt_pdata(2:end);
pind = find(pddfilt_pdata <= 0);
csind = [1; pind(1:end)+1];     % long column of cast start time
ceind = [pind(1:end); length(pdata)];   % cast endings
thind = [];
for i = 1:length(csind)
    if ceind(i)-csind(i) < 1000  % not a profile
        thind = [thind i]; 
    end
end
csind(thind) = [];
ceind(thind) = [];
fprintf('identify %i profile \n',round(length(csind)/2))


% we start at the top when the fallrate is higher than the criterium
% and look for the next time the speed is lower than that criteria to end
% the cast.
% we iterate the process to define all the casts
Start_ind    =  find(pdata>min_depth,1,'first');
nb_down   =  1;
nb_up     =  1;
do_it        =  0;

while (do_it==0)
    End_ind=Start_ind+find(pdata(Start_ind+1:end)<min_depth,1,'first');
    if ~isempty(End_ind)
        down{nb_down}=Start_ind+(1:find(pdata(Start_ind:End_ind)==max(pdata(Start_ind:End_ind))));
        nb_down=nb_down+1;
        up{nb_up}=Start_ind+find(pdata(Start_ind:End_ind)==max(pdata(Start_ind:End_ind)))+1:End_ind;
        nb_up=nb_up+1;
        Start_ind =  End_ind+1+find(pdata(End_ind+1:end)>min_depth,1,'first');
    else
        do_it=1;
    end
    if mod(nb_up,10)==0
        fprintf('Upcast %i\n',nb_up)
        fprintf('Downcast %i \n',nb_up)
    end
end


%once we have the index defining the casts we split the data
dataup=cellfun(@(x) structfun(@(y) y(x),data,'un',0),up,'un',0);
datadown=cellfun(@(x) structfun(@(y) y(x),data,'un',0),down,'un',0);


% select only cast with more than 10 points. 10 points is arbitrary
indPup=find(cellfun(@(x) x.P(1)-x.P(end),dataup)>crit_depth);


indPdown=find(cellfun(@(x) x.P(end)-x.P(1),datadown)>crit_depth);

datadown=datadown(indPdown);
dataup=dataup(indPup);
down=down(indPdown);
up=up(indPup);

% put info back
if isfield(data,'info')
    for i = 1:length(dataup)
        datadown{i}.info=info;
        dataup{i}.info=info;
    end
end
clear i

if isfield(data,'sn')
for i = 1:length(dataup)
        datadown{i}.sn=sn;
        dataup{i}.sn=sn;
end
end
clear i

% plot pressure and highlight up /down casts in red/green
figure
% plot(tdata,pdata)
hold on
for i=1:length(up)
    if isfield(data,'ctdtime')
        plot(dataup{i}.ctdtime,dataup{i}.P,'r')
    else
        plot(dataup{i}.time,dataup{i}.P,'r')
    end
end
for i=1:length(down)
    if isfield(data,'ctdtime')
       plot(datadown{i}.ctdtime,datadown{i}.P,'b')
    else
       plot(datadown{i}.time,datadown{i}.P,'b')
    end
end
axis ij
% 
% %MHA: plot ocean style
% axis ij

% print('-dpng2',[CTD.CTDpath 'Profiles_Pr.png'])


% % do we wnat to save or change the speed and filter criteria
% answer1=input('save? (yes,no)','s');
% 
% %save
% switch answer1
%     case 'yes'
%         
%         CTDProfiles.up=up;
%         CTDProfiles.down=down;
%         CTDProfiles.dataup=dataup;
%         CTDProfiles.datadown=datadown;
%         CTDProfiles.info=info;
%         filepath=fullfile(CTD.L1path,['Profiles_' CTD.deployment '.mat']);
%         fprintf('Saving data in %s \n',filepath)
%         save(filepath,'CTDProfiles','-v7.3');
%         
%         CTD.nbprofileup=numel(CTDProfiles.up);
%         CTD.nbprofiledown=numel(CTDProfiles.down);
%         CTD.maxdepth=max(cellfun(@(x) max(x.P),CTDProfiles.dataup));
%         save(fullfile(CTD.L1path,'Meta_Data.mat'),'Meta_Data')
% end
% 


end