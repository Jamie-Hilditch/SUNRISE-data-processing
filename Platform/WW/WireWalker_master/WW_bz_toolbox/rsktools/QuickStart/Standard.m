%% RSKtools for Matlab access to RBR data
% RSKtools v2.2.0;
% RBR Ltd. Ottawa ON, Canada;
% support@rbr-global.com;
% 2018-01-25

%% Introduction 
% |RSKtools| provides convenient functions for data extraction (e.g.,
% extracting profiles from a continuous dataset) and visualisation
% (e.g., plotting individual profiles). As of v2.0.0, RSKtools
% includes a suite of functions to perform routine processing steps to
% enhance the data quality (see the Resources section for information
% on RSKtools post-processing functions). We are continually expanding
% RSKtools, and please feel free to suggest improvements and new
% features.

%% Installing
% The latest stable version of |RSKtools| can be found at
% <http://www.rbr-global.com/support/matlab-tools>.
% 
% * Download and unzip the archive (to |~/matlab/RSKtools|, for instance) 
% * Add the folder to your path from the command line (|addpath ~/matlab/RSKtools|) or launch the path editor gui (|pathtool|). 
% * type |help RSKtools| to get an overview and take a look at the examples.

  
%% Examples of use
% <html><h3>Loading files</h3></html>
% 
% The first step is to make a connection to the RSK file with
% |RSKopen|. Note that |RSKopen| does not actually read the data;
% instead it reads a "thumbnail" of the data, which is up to 4000
% samples long. The structure returned after opening an RSK looks
% something like:

file = '../sample.rsk';
rsk = RSKopen(file)

%%
% To read the full dataset, use the |RSKreaddata| function.
% RSKreaddata will read the full dataset by default.  Because RSK
% files can store a large amount of data, it may be preferable to read
% a subset of the data, specified using a start and end time (in
% Matlab |datenum| format, which is defined as the number of days
% since January 0, 0000).

t1 = datenum(2014, 05, 03);
t2 = datenum(2014, 05, 04);
rsk = RSKreaddata(rsk, 't1', t1, 't2', t2);

%%
% Note that the logger data can be found in the structure at:

rsk.data        

%%
% where |rsk.data.tstamp| contains the sample timestamps in Matlab
% datenum format, and |rsk.data.values| contains the data.  Each
% column in |rsk.data.values| contains data from a different channel.
% The channel names and units for each column in |data| are:

[{rsk.channels.longName}' {rsk.channels.units}']

%% 
% In this particular example, Practical Salinity can be derived from
% conductivity, temperature, and pressure because the file comes from
% a "CTD"-type instrument.  |RSKderivesalinity| is a wrapper for the
% TEOS-10 GSW function |gsw_SP_from_C|, and it adds a new channel
% called |Salinity| as a column in |rsk.data.values|.  The TEOS-10 GSW
% Matlab toolbox is freely available from
% <http://teos-10.org/software.htm>.  It is good practice to derive
% sea pressure first in case you want to specify a custom value of
% atmospheric pressure, otherwise the nominal value of 10.1325 dbar is
% used.

rsk = RSKderiveseapressure(rsk);
rsk = RSKderivesalinity(rsk);


%% Working with profiles
% Profiling loggers with recent versions of firmware can detect and
% record profile upcast and downcast "events" automatically. The
% function |RSKreadprofiles| uses the profile event time stamps to
% organize the data into profiles. Then, a plot of the profiles can be
% made very easily using the |RSKplotprofiles| function.
%
% If profiles have not been detected by the logger or Ruskin, or if
% the profile time stamps do not correctly parse the data into
% profiles, the function |RSKfindprofiles| can be used. The
% |pressureThreshold| argument, which determines the pressure reversal
% required to trigger a new profile, and the |conductivityThreshold|
% argument, which determines if the logger is out of the water, can be
% adjusted to improve profile detection when the profiles were very
% shallow, or if the water was very fresh.
%
% Note: Salinity, sea pressure, and other channels derived by RSKtools
% should be derived again after using |RSKreadprofiles| because in
% some circumstances it will go back and read raw data from the RSK
% file.

% load the upcast and downcast of profiles 6 to 8
rsk = RSKreadprofiles(rsk, 'profile', 6:8, 'direction', 'both');
rsk = RSKderiveseapressure(rsk);
rsk = RSKderivesalinity(rsk);

% plot the upcasts of temperature, salinity, and chlorophyll
handles = RSKplotprofiles(rsk, 'channel', {'temperature','salinity','chlorophyll'}, 'direction', 'up');


%% Customising plots
% The plotting functions return a handle enabling access to the line
% objects. The output is a matrix containing a column for each channel
% subplot and a row for each profile.
handles

% To increase the line width of the first profile in all subplots
set(handles(1,:),{'linewidth'},{3});


%% Accessing individual channels and profiles
% The channel data is stored in |rsk.data|.  If the data was parsed
% into profiles, |data| is a 1xN structure array, where each element
% is an upcast or downcast from a single profile containing an array
% of time stamps and a matrix of channel data.  |RSKtools| has
% functions to access the data from particular channels and profiles.
% For example, to access the time stamps, sea pressure, temperature,
% and dissolved O2 from the upcast of the 1st profile, run:

profind = getdataindex(rsk,'direction','up','profile',1);
tempcol = getchannelindex(rsk,'temperature');
o2col   = getchannelindex(rsk,'dissolved o2');
prescol = getchannelindex(rsk,'sea pressure');

time        = rsk.data(profind).tstamp;
seapressure = rsk.data(profind).values(:,prescol);
temperature = rsk.data(profind).values(:,tempcol);
o2          = rsk.data(profind).values(:,o2col);


%% Other Resources
% We recommend reading:
%
% * The <https://docs.rbr-global.com/rsktools RSKtools on-line user
% manual> for detailed RSKtools function documentation.
%
% * The <http://rbr-global.com/wp-content/uploads/2017/08/VignettePostProcessing.pdf
% RSKtools post-processing guide> for an introduction on how to
% process RBR profiles with RSKtools.  The post-processing suite
% contains, among other things, functions to low-pass filter, align,
% de-spike data, trim the data, and write CSV files.


%% About this document
% This document was created using
% <http://www.mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html
% Matlab(TM) Markup Publishing>. To publish it as an HTML page, run the
% command:
%%
% 
%   publish('Standard.m');

%%
% See |help publish| for more document export options.