function [pos,dist] = latlon2xy(lat,lon,lat_c,lon_c)

radius=6373.19*1e3;

dy = radius*pi/180*1; % 1 degree
dx = radius*cosd(lat_c)*pi/180*1;

pos =[(lon-lon_c)*dx,(lat-lat_c)*dy]';

dist = vecnorm(diff(pos'),2,2);

return