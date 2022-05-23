function data = parse_rbr_concerto(f_in)
    data_path

    raw_mat = RSKopen(f_in);
    raw_mat = RSKreaddata(raw_mat);
    
    [raw_mat,~] = RSKcorrecthold(raw_mat);
    [~,temp1] = min(abs(raw_mat.data.tstamp - raw_mat.data.tstamp(1) - 1/1440));
    raw_mat = RSKderiveseapressure(raw_mat,'patm',nanmean(raw_mat.data.values(1:temp1,getchannelindex(raw_mat, 'Pressure'))));
    raw_mat = RSKderivedepth(raw_mat,'latitude',lat_c);
    raw_mat = RSKderivevelocity(raw_mat);
    
    CT_lag = RSKcalculateCTlag(raw_mat);
    raw_mat = RSKalignchannel(raw_mat, 'channel', 'Conductivity', 'lag', CT_lag);
    
    raw_mat = RSKderivesalinity(raw_mat);
    raw_mat = RSKderivetheta(raw_mat,'latitude',lat_c,'longitude',lon_c);
    raw_mat = RSKderiveSA(raw_mat,'latitude',lat_c,'longitude',lon_c);
    raw_mat = RSKderivesigma(raw_mat,'latitude',lat_c,'longitude',lon_c);

    data = struct();

    CTD_name = {'Pressure','Depth','Conductivity','Salinity','Absolute Salinity','Temperature','Potential Temperature',...
    'Density Anomaly'};
    CTD_short = {'P','depth','C','SP','SA','T','theta','sigma'};

    for k = 1 :length(CTD_name)
        idx = find(ismember({raw_mat.channels.longName},CTD_name{k}));
        data.(CTD_short{k}) = raw_mat.data.values(:,idx);
    end

    data.dn = raw_mat.data.tstamp;
end
