function RBRgrid = mod_telemetry_gridctd(dataup,dz,zmax)

% load(fullfile(WWmeta.mappath,['Profiles_' WWmeta.name_rbr]),'RBRprofiles')
for i = 1:length(dataup)
if isfield(dataup{i},'info')
    info=dataup{i}.info;
    dataup{i}=rmfield(dataup{i},'info');
end

if isfield(dataup{i},'sn')
    sn=dataup{i}.sn;
    dataup{i}=rmfield(dataup{i},'sn');
end

if isfield(dataup{i},'download')
    download=dataup{i}.download;
    dataup{i}=rmfield(dataup{i},'download');
end
end

%get the normal upcast (mean P of the upcast ~ median P of all the mean P)

RBRprofiles = dataup;
Prbr=cellfun(@(x) mean(x.P),RBRprofiles);
critp= max(cellfun(@(x) max(x.P),RBRprofiles))-.5*std(Prbr);
critm= min(cellfun(@(x) min(x.P),RBRprofiles))+.5*std(Prbr);


timerbr=cellfun(@(x) mean(x.time),RBRprofiles);
timerbrOK=timerbr(Prbr>critm & Prbr<critp);
indOK=(Prbr>critm & Prbr<critp);
RBRprofiles=RBRprofiles(indOK);

if nargin==2
    zaxis=0:dz:max(cellfun(@(x) max(x.P),RBRprofiles));
elseif nargin==3
    zaxis=0:.25:100';
else
    zaxis=0:.25:max(cellfun(@(x) max(x.P),RBRprofiles));
end

Z=length(zaxis);

fields=fieldnames(RBRprofiles{1});
for f=1:length(fields)
    wh_field=fields{f};
    if ~strcmp(wh_field,'info')
        RBRgrid.(wh_field)=zeros([Z,sum(indOK)]);
        for t=1:length(timerbrOK)
            F=RBRprofiles{t}.(wh_field);
            [Psort,I]=sort(RBRprofiles{t}.P,'descend');
            P_temp=(interp_sort(Psort));
            RBRgrid.(wh_field)(:,t)=interp1(P_temp,F(I),zaxis);
        end
    end
end
RBRgrid.z=zaxis'*ones(1,size(RBRgrid.T,2));
RBRgrid.time=ones(size(RBRgrid.T,1),1)*timerbrOK/86400/1000+719529;   % in days
RBRgrid.info = info;
RBRgrid.sn = sn;
% RBRgrid.download = download;




% %% clean n2 field
% n2=RBRgrid.n2;
% n2(n2<=0)=nan;
% [TT,ZZ]=meshgrid(RBRgrid.time,RBRgrid.z);
% %TT=TT(:);ZZ=ZZ(:);
% F=scatteredInterpolant(TT(~isnan(n2)),ZZ(~isnan(n2)),n2(~isnan(n2)));
% nn2=F(TT(:),ZZ(:));
% nn2=reshape(nn2,[length(RBRgrid.z),length(RBRgrid.time)]);
% nn2(nn2<=0)=1e-10;
% RBRgrid.n2=nn2;
% 
% RBRgrid.z=zaxis;
% RBRgrid.time=timerbrOK;
% RBRgrid.info=RBRprofiles{1}.info;

% save(fullfile(WWmeta.mappath,[WWmeta.WW_name '_grid.mat']),'RBRgrid')
% save(fullfile(WWmeta.mappath,['Grid_' WWmeta.name_rbr]),'RBRgrid')

%% TODO: uncomment and improve when Nortek processing is reading
% if exist(fullfile([WWmeta.WWpath WWmeta.WW_name '_grid.mat']),'file')
%     load(fullfile([WWmeta.WWpath WWmeta.WW_name '_grid.mat']),'AQDgrid')
%     save(fullfile([WWmeta.WWpath WWmeta.WW_name '_grid.mat']),'RBRgrid','AQDgrid')
% else
%     save(fullfile([WWmeta.WWpath WWmeta.WW_name '_grid.mat']),'RBRgrid')
% end

