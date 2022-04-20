function [RBRprofiles_up, RBRprofiles_down]=create_profiles_rbr(CTD)
% identify upcast and downcast and then save
%
% name_file = [WWmeta.name_aqd,'_',num2str(q),'_',num2str(q+num-1)];
% load([WWmeta.matpath name_file '.mat'])  % load data
% eval(['[up,down,dataup,datadown] = get_aqd_2G(CTD);']) % identify upcast and downcast
[up,down,dataup,datadown] = get_rbr_2G(CTD);
% if isfield(WWmeta,'figpath')  % check whehter need to save the figures
%     print([WWmeta.figpath,name_file,'.png'],'-r300', '-dpng');  % save the figure (optional)
% end

% generate upcast
dup=diff(up);
ind_prof=find(dup>1);
ind_prof = [1;ind_prof;length(up)]; % change by BZheng on Feb 15, 2020, orgignal: without this command
RBRprofiles_up=struct([]);
fields=fieldnames(dataup);
tdata=dataup.time;
for i=1:length(ind_prof)-1   
    for f=1:length(fields)
        wh_field=fields{f};
        if (length(tdata)==length(dataup.(wh_field)))
            switch length(size(dataup.(wh_field)))
                case 2
                   RBRprofiles_up{i}.(wh_field)=dataup.(wh_field)(ind_prof(i)+1:ind_prof(i+1),:);
                case 3
                   RBRprofiles_up{i}.(wh_field)=dataup.(wh_field)(ind_prof(i)+1:ind_prof(i+1),:,:);
                case 4
                   RBRprofiles_up{i}.(wh_field)=dataup.(wh_field)(ind_prof(i)+1:ind_prof(i+1),:,:,:);
            end
        end
    end
end

% generate downcast
ddown=diff(down);
ind_prof=find(ddown>1);  % find gap, which is an indication of profile
ind_prof = [1;ind_prof;length(down)]; % change by BZheng on Feb 15, 2020, orgignal: without this command
RBRprofiles_down=struct([]);
fields=fieldnames(datadown);
tdata=datadown.time;
for i=1:length(ind_prof)-1   
    for f=1:length(fields)
        wh_field=fields{f};
        if (length(tdata)==length(datadown.(wh_field)))
            switch length(size(datadown.(wh_field)))
                case 2
                   RBRprofiles_down{i}.(wh_field)=datadown.(wh_field)(ind_prof(i)+1:ind_prof(i+1),:);
                case 3
                   RBRprofiles_down{i}.(wh_field)=datadown.(wh_field)(ind_prof(i)+1:ind_prof(i+1),:,:);
                case 4
                   RBRprofiles_down{i}.(wh_field)=datadown.(wh_field)(ind_prof(i)+1:ind_prof(i+1),:,:,:);
            end
        end
    end
end

% save([WWmeta.propath 'Profiles_upcast_' name_file],'AQDprofiles_up')
% save([WWmeta.propath 'Profiles_downcast_' name_file],'AQDprofiles_down')


