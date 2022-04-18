% Program written to incorporate CF warm layer code in coare35vn bulk
% algorithm
% Ludovic Bariteau & Jim Edson
% 1st version: March 2013

%clear all;
%close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grab the 10 minute file from DYNAMO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd c:\data\cwf\dropbox\matlabstf\cwf\bulkalg\cor3_5;
load Revelle10minutesLeg3_r3; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs Required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Yday           %decimal yearday
%Ur             %water-relative wind vector magnitude, m/s
%Tsea           %Sea temperature, C at ts_depth 
%Tair           %air temperature, C
%Pair           %air temperature, C
%RH             %relative humidity (%)
%Lat            %latitude, deg
%Lon;           %longitude, deg
%Solar          %downward solar flux, W/m^2
%IR             %downward IR flux, W/m^2
%Rainrate       %rainrate, mm/hr
%Ss             %sea surface salinity (PSU), use 35 if unknown
%cp             %phase speed of dominant waves (m/s), set to NaN if unknown  
%sigH           %significant wave height (m, set of NaN if unknown
%zu, zt, zq     %heights of the observations (m)
%zrf_u, zrf_t, zrf_q  reference height for profile.  Use this to compare observations at different heights  
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zu=10*ones(1,length(U10));;             
zt=10*ones(1,length(U10));;             
zq=10*ones(1,length(U10));;             
Yday=yday;          %decimal yearday Jan 1 at noon = 1.5
Ur=Ur10;            %relative wind speed
U=U10;              %true wind speed
Pair=Pair10;        %pressure, mb
Tair=T10;           %air temperature, C
RH=RH10;            %relative humidity, %
Tsea=Tsea;          %bulk sea temperature, C
ts_depth=2;         %bulk sea temperature sensor depth
Solar=-Solardn;     %downward solar flux, W/m^2
IR=-IRdn;           %downward IR flux, W/m^2
Rainrate=P;         %rainrate, mm/hr
Lat=Lat;            %latitude, deg
Lon=Lon;            %longitude, deg
zi=600;%*ones(1,length(U));             %inversion ht
Ss=Sal5;             %surface salinity, ppt
zrf_u=10*ones(1,length(U));zrf_t=10*ones(1,length(U));zrf_q=10*ones(1,length(U));%reference heights for means, m
cp=NaN*U;            %wave peak phase speed, m/s
sigH=NaN*U;          %wave significant height, m



%%%%%%%%%%%%%   vectorized coare algorithm fluxes computed with TSG at
%%%%%%%%%%%%%   depth of ts_depth
A36=coare36vn_zrf(U,zu,Tair,zt,RH,zq,Pair,Tsea5,Solar,IR,Lat,zi,Rainrate,Ss,cp,sigH,zrf_u,zrf_t,zrf_q);
figure;plot(yday,A36(:,4),yday,A36(:,3));xlabel('Yearday (2011)');ylabel('Heat Flux (W/m^2)');title('Bulk heat fluxes from DYNAMO');legend('Hl','Hs');
%A35=coare35vn(U,zu,Tair,zt,RH,zq,Pair,Tsea5,Solar,IR,Lat,zi,Rainrate,cp,sigH);%% version 35

%************************************************************
%  This test program uses the TSG temperature at
%  approximately 2m to compute the depth and amplitude of 
%  the warm layer.  The heating at this depth is added to 
%  TSG and the cool-skin affect is removed in the 
%  coare36vn algorithm.  Therefore, fluxes in the returned 
%  B array have been computed based on a TSG measurements 
%  that have been corrected for the warm layer and cool skin.
%  The sea snake temperature is consider to be the 'true' surface temperature 
%************************************************************

B=coare36vnWarm(yday,U,zu,Tair,zt,RH,zq,Pair,Tsea5,Solar,IR,Lat,Lon,zi,Rainrate,ts_depth,Ss,cp,sigH,zrf_u,zrf_t,zrf_q);

%Outputs
dt_wrmx=B(:,44)'; %warming across entire warm layer deg.C
tk_pwpx=B(:,45)';%warm layer thickness m
dseax=B(:,46)';%heating at selected depth

%************************************************************
%  The following shows how to use the warm layer variables to
%  compare the temperature measured at another depth, i.e.,
%  the TSG. 
%  If the sensor is beneath the warm layer depth, the total 
%  temperature change (the warm layer amplitude) is added to 
%  the sensor. 
%  Note that these are the values just below the surface 
%  without the cool skin correction.
%************************************************************

tsg_depthx=2;                       %Approximate depth of TSG
dseagx=dt_wrmx.*tsg_depthx./tk_pwpx;%Linear heating profile is assumed
i=find(tk_pwpx<tsg_depthx);         %Use the total warming if the sensor is below the warm layer.
dseagx(i)=dt_wrmx(i);               %
tsgx=Tsea5+dseagx;                  %Add warm layer correction to adjust temps to just below surface values.
tsx=Tsea;                           %Assume sea snake is true SST
dtwgx=dt_wrmx-dseagx;

figure(1);
plot(yday,Tsea,'r',yday,Tsea5,'g','linewidth',2);xlabel('Yearday (2011)');ylabel('Temperature (C)');
hold
plot(yday,tsgx,'o')
axis([315 325 29 33]);
legend('Sea Snake','ThermoSalinoGraph','TSG Corrected','Location','Northeast')
hold off
