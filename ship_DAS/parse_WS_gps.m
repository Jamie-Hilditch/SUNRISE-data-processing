%% parse_WS_gps.m
% Usage: parse_WS_gps
% Description: Read navigational data from files in the Walton Smith's 'POSMV'
%              directory and create a single output .mat file.
% Inputs: None
% Outputs: None
%
% Author: Dylan Winters
% Created: 2022-01-06
clear all, close all

dir_in = fullfile('../../data/raw/ship/WaltonSmith/POSMV')
tbl = parse_WS_DAS(dir_in,'GGA');

% Create output structure
out.dn = datenum(tbl.ComputerDate + tbl.ComputerTime); % date & time to datenum
out.lat = tbl.LatDeg + 1/60*tbl.LatMin;                % deg/min to decimal deg
out.lon = -(tbl.LonDeg + 1/60*tbl.LonMin);             % deg/min W to decimal deg E

% Clear the large table
clear tbl

%% Heading: in files containing HDT messages
tbl = parse_WS_DAS(dir_in,'HDT');

% Add to output structure
dn = datenum(tbl.ComputerDate + tbl.ComputerTime); % date & time to datenum
h = tbl.TrueHeading;

% Interpolate heading onto lat/lon time vector
hi = exp(1i*pi/180*h); % make imaginary
[~,iu] = unique(dn);
hi = interp1(dn(iu),hi(iu),out.dn,'spline'); % interpolate
out.heading = mod(180/pi*angle(hi),360);     % convert back to degrees

save('WS_nav.mat','-struct','out');
