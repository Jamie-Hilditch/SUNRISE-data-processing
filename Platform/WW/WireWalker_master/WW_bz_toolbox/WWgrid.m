function WWgrid(WWmeta,zgrid)
% grid ww profile data
% bz, june 15, 2021

load([WWmeta.propath,'CTDprofiles.mat']);

%get the normal upcast (mean P of the upcast ~ median P of all the mean P)

Prbr=cellfun(@(x) mean(x.P),RBRprofiles);
critp= max(cellfun(@(x) max(x.P),RBRprofiles))-.5*std(Prbr);
critm= min(cellfun(@(x) min(x.P),RBRprofiles))+.5*std(Prbr);

timerbr=cellfun(@(x) mean(x.time),RBRprofiles);
timerbrOK=timerbr(Prbr>critm & Prbr<critp);
indOK=(Prbr>critm & Prbr<critp);
RBRprofiles=RBRprofiles(indOK);

% combine all upcasts
fields=fieldnames(RBRprofiles{1});
idx = size(length(RBRprofiles),2);
for f=1:length(fields)
    wh_field=fields{f};
    if ~strcmp(wh_field,'info')
        RBRgrid.(wh_field)=[];
        for t=1:length(RBRprofiles)
            F=RBRprofiles{t}.(wh_field);
            RBRgrid.(wh_field)=[RBRgrid.(wh_field); F];
            if f == length(fields)
                idx(t,2) = length(RBRgrid.(wh_field));
                idx(t,1) = length(RBRgrid.(wh_field)) - length(F) + 1;
            end
            
        end
    end
    
end
RBRgrid.idx = idx;

%% put gridded product into std_profiles struct
if nargin==2
    zaxis=0:zgrid:max(cellfun(@(x) max(x.P),RBRprofiles));
else
    zaxis=0:.25:max(cellfun(@(x) max(x.P),RBRprofiles));
end

Z=length(zaxis);

fields=fieldnames(RBRprofiles{1});
for f=1:length(fields)
    wh_field=fields{f};
    if ~strcmp(wh_field,'info')
        RBRgrid.std_profiles.(wh_field)=zeros([Z,sum(indOK)]);
        for t=1:length(timerbrOK)
            F=RBRprofiles{t}.(wh_field);
            [Psort,I]=sort(RBRprofiles{t}.P,'descend');
            P_temp=(interp_sort(Psort));
            RBRgrid.std_profiles.(wh_field)(:,t)=interp1(P_temp,F(I),zaxis);
        end
    end
end
RBRgrid.std_profiles.z=zaxis';  % column array
RBRgrid.std_profiles.time=timerbrOK;
% RBRgrid.std_profiles.info=RBRprofiles{1}.info;

%%
save([WWmeta.gridpath,'CTDgrid.mat'],'RBRgrid');
end