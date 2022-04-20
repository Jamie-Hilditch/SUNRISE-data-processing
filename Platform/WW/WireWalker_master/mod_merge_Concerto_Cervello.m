function CTD = mod_merge_Concerto_Cervello(Cervello,Concerto)
% function to interpolate Cervello timestamps on to Concerto timestamps

CTD = Concerto;

% interpolate
CTD.lat = interp1(Cervello.time,Cervello.lat,Concerto.time);
CTD.lon = interp1(Cervello.time,Cervello.lon,Concerto.time);

end