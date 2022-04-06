%% load_uhdas.m
% Usage: adcp = load_uhdas(dir_in)
% Description: Read ADCP data processed by UHDAS.
% Inputs: dir_in - string, path to UHDAS directory
% Outputs: adcp - struct, looks something like this:
%            dn: [1x25230 double]
%           lat: [1x25230 double]
%           lon: [1x25230 double]
%       heading: [1x25230 double]
%         depth: [40x1 double]
%             z: [40x1 double]
%             u: [40x25230 double]
%        u_ship: [1x25230 double]
%             v: [40x25230 double]
%        v_ship: [1x25230 double]
%             w: [40x25230 double]
% Author: Dylan Winters
% Created: 2021-12-15

% NOTE: These are not the final variable names

function adcp = load_uhdas(dir_in)

% Directory containing matfiles
mdir = fullfile(dir_in,'contour');

% Initialize output structure
adcp = struct();

% Load time & position information
f_in = fullfile(mdir,'allbins_other.mat');
tmp = load(f_in);
adcp.dn = datenum(tmp.TIME)';
adcp.lat = tmp.LAT_END';
adcp.lon = tmp.LON_END';
nb = tmp.NBINS;
adcp.heading = tmp.HEADING';

% Load depth information
f_in = fullfile(mdir,'contour_xy.mat');
tmp = load(f_in);
adcp.depth = tmp.z(1:nb);
adcp.z = - adcp.depth;

% Load velocity
vels = {'u','v','w'};
for i = 1:length(vels)
    f_in = fullfile(mdir,sprintf('allbins_%s.mat',vels{i}));
    tmp = load(f_in);
    adcp.(vels{i}) = tmp.(upper(vels{i}));

    % Remove ship velocity from measured velocity
    vel_ship = [upper(vels{i}) '_SHIP'];
    if isfield(tmp,vel_ship);
        adcp.(lower(vel_ship)) = tmp.(vel_ship);
        adcp.(vels{i}) = adcp.(vels{i}) + tmp.(vel_ship);
    end
end
