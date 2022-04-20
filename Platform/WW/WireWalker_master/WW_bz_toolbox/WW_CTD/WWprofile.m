function WWprofile(WWmeta,theshold)
% separate raw data into profiles
% bz, june 15, 2021
load([WWmeta.matpath,'CTDall.mat']);
[up,down,dataup,datadown] = get_ctd_2G(CTDall,theshold);  % find upcast/downcast, pay attention to this threshold value

dup=diff(up);
ind_prof=find(dup>1);
ind_prof = [0; ind_prof; length(up)];  % add the first index and the last index
RBRprofiles=struct([]);
fields=fieldnames(dataup);
tdata=dataup.time;
for i=1:length(ind_prof)-1
    for f=1:length(fields)
        wh_field=fields{f};
        if (length(tdata)==length(dataup.(wh_field)))
            RBRprofiles{i}.(wh_field)=dataup.(wh_field)(ind_prof(i)+1:ind_prof(i+1));
%             RBRprofiles{i}.info=CTDall.info;
        end
    end
end

index=find(abs(cellfun(@(x) x.P(end)-x.P(1),RBRprofiles))<5);
RBRprofiles(index) = [];

%%
save([WWmeta.propath,'CTDprofiles.mat'],'RBRprofiles');
end