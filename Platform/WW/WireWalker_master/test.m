for q = 1:variables.NUM_combining_files:442
    if q+variables.NUM_combining_files-1>442 % if length of the index is longer than the size of the raw data, then the ending index is set to be the lenght of the raw data file
        num = 442-q+1;
    else
        num = variables.NUM_combining_files;
    end
    
%     merge_signature(WWmeta,q,num);    % merge separate .mat files from ADCP output and then form a group
    create_profiles(WWmeta,q,num,200);  % will chunk ww profiles into upcast and downcast, and save them separately
    disp(['current file location: ',num2str(q),'_',num2str(q+num-1)])  % show where we are
end
disp('identify profiles: finished')