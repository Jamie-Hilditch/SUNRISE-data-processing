function B=coare36vnWarm(yday,Ur,zu,Tair,zt,RH,zq,Pair,Tsea,Solar,IR,Lat,Lon,zi,Rainrate,ts_depth,Ss,cp,sigH,zrf_u,zrf_t,zrf_q);

disp('WarmCoolLayer')

%***********   input data **************
%   yday=       day-of-year
%	Ur=			wind vector maginitude (m/s) relative to water at height zu
%	zu=			height (m) of wind measurement
%	Tair=		air temp (degC)at height zt
%	zt=			height (m) of air temperature measurement
%	RH=			relative humidity (%) at height zq
%	zq=			height (m) of air humidity measurement
%	Pair=		air pressure (mb)
%	Tsea=	    bulk surface sea temp (degC) at ts_depth
%	Solar=		downward solar flux (w/m^2) defined positive down
%	IR=			downward IR flux (w/m^2) defined positive down
%	Lat=		latitude (deg N=+)
%	Lon=		longitude (deg E=+)
%   zi=         inversion height (m)
%	Rainrate=	rain rate (mm/hr)
%	ts_depth	depth (m) of water temperature measurement
%   Ss =        sea surface salinity (PSU)
%   cp =        phase speed of dominant waves (m/s)  
%   sigH =      significant wave height (m)
%   zu, zt, zq heights of the observations (m)
%   zrf_u, zrf_t, zrf_q  reference height for profile.  Use this to compare observations at different heights  
%
%  
%********** output data  ***************
%Outputs
%From coare36vn_zrf
%Output:   A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq  Cd  Ch  Ce  L zet dterx dqerx tkt Urf Trf Qrf RHrf UrfN Rnl  Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis];
%              1   2   3   4   5    6     7    8   9  10  11   12  13  14  15 16  17   18   19    20  21  22  23  24   25   26  27  28  29  30  31    32     33     34   35 36  37   38  39  40  41  42  43 
%where
%
%   usr = friction velocity that includes gustiness (m/s)
%   tau = wind stress (N/m^2)
%   hsb = sensible heat flux into ocean (W/m^2)
%   hlb = latent heat flux into ocean (W/m^2)
%   hbb = buoyany flux into ocean (W/m^2)
%   hsbb = "sonic" buoyancy flux measured directly by sonic anemometer 
%   tsr = temperature scaling parameter (K)
%   qsr = specific humidity scaling parameter (g/Kg)
%   zot = thermal roughness length (m)
%   zoq = moisture roughness length (m)
%   Cd = wind stress transfer (drag) coefficient at height zu   
%   Ch = sensible heat transfer coefficient (Stanton number) at height zu   
%   Ce = latent heat transfer coefficient (Dalton number) at height zu
%    L = Obukhov length scale (m) 
%  zet = Monin-Obukhov stability parameter zu/L 
% dter = cool-skin temperature depression (degC)
% dqer = cool-skin humidity depression (degC)
%  tkt = cool-skin thickness (m)
%  Urf = wind speed at reference height (user can select height below)
%  Tfr = temperature at reference height
%  Qfr = specific humidity at reference height
% RHfr = relative humidity at reference height
% UrfN = neutral value of wind speed at reference height
%  Rnl = Upwelling IR radiation computed by COARE
%   Le = latent heat of vaporization
% rhoa = density of air
%   UN = neutral value of wind speed at zu
%  U10 = wind speed adjusted to 10 m
% UN10 = neutral value of wind speed at 10m
%Cdn_10 = neutral value of drag coefficient at 10m    
%Chn_10 = neutral value of Stanton number at 10m    
%Cen_10 = neutral value of Dalton number at 10m    
%Rf     = Rain heat flux (W/m^2)
%Qs     = surface specific humidity (g/g)
%Evap   = evaporation rate (mm/h)
%T10    = air temperature at 10m
%Q10    = air specific humidity at 10m
%RH10   = air relative humidity at 10m
%uq     = gustiness velocity (m/s)
%Whf    = whitecap fraction
%Edis   = energy dissipated by wave breaking (W/m^2)



%From WarmLayer
% 44: dt_wrm - warming across entire warm layer degC
% 45: tk_pwp - warm layer thickness m
% 46: dsea   - dT due to warming at depth of Tsea such that Tsea_true = Tsea + dsea

jcool=1;  %0=no cool skin calc, 1=do cool skin calc
icount=1;
%*********************  housekeep variables  ********
% Call coare35vn to get initial flux values
Bx=coare36vn_zrf(Ur(1),zu(1),Tair(1),zt(1),RH(1),zq(1),Pair(1),Tsea(1),Solar(1),IR(1),Lat(1),zi(1),Rainrate(1),Ss(1),cp(1),sigH(1),zrf_u(1),zrf_t(1),zrf_q(1));
tau_old=Bx(2);  %stress
hs_old=Bx(3);   %sensible heat flux
hl_old=Bx(7);   %latent heat flux - Webb corrected
dter=Bx(17);    %cool skin
RF_old=Bx(34);  %rain heat flux

qcol_ac=0;      %accumulates heat from integral
tau_ac=0;       %accumulates stress from integral 
dt_wrm=0;       %total warming (amplitude) in warm layer
max_pwp=19;     %maximum depth of warm layer (adjustable)
tk_pwp=max_pwp; %initial depth set to max value
dsea=0;         %dT initially set to 0
q_pwp=0;        %total heat absorped in warm layer
fxp=.5;         %initial value of solar flux absorption

rich=.65;       %critical Richardson number

jtime=0;
jamset=0;
jump=1;

%*******************  set constants  ****************
tdk=273.16;   %Converts to Kelvin
Rgas=287.1;   %Universal gas constant
cpa=1004.67;  %Specific heat of air at constant pressure
cpw=4000;     %Specific heat of water
rhow=1022;    %Density of water
visw=1e-6;    %Viscosity of water

be=0.026;
tcw=0.6;

%**********************************************************
%******************  setup read data loop  ****************

[Press,Tseak,Tairk,Qsatsea,Qsat,Qair,Rhoair,Rhodry]=scalarv(Pair,Tsea,Tair,RH);

nx=length(yday);        %# of lines of data

for ibg = 1:nx 			%major read loop
    yd=yday(ibg);       %yearday
    P=Pair(ibg);        %air pressure
    u=Ur(ibg);          %wind speed
    tsea=Tsea(ibg);     %bulk sea surface temp
    t=Tair(ibg);        %air temp
    qs=Qsatsea(ibg);    %bulk sea surface humidity
    q=Qair(ibg);        %specific humidity  
    rh=RH(ibg);         %relative humidity
    Rs=Solar(ibg);      %downward solar flux (positive down)
    Rl=IR(ibg);         %doward IR flux (positive down)
    rain=Rainrate(ibg); %rain rate
    grav=grv(Lat(ibg)); %gravity
    latx=Lat(ibg);      %latitude
    lonx=Lon(ibg);      %longitude
    rhoa=Rhoair(ibg);   %air density
    
    %*****  variables for warm layer  ***
    Rnl=.97*(5.67e-8*(tsea-dter*jcool+tdk)^4-Rl); %Net IR
    Rns=.945*Rs;                                  %Net Solar
    cpv=cpa*(1+0.84*q/1000);
    visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t);
    Al=2.1e-5*(tsea+3.2)^0.79;
    ctd1=sqrt(2*rich*cpw/(Al*grav*rhow));       %mess-o-constants 1
    ctd2=sqrt(2*Al*grav/(rich*rhow))/(cpw^1.5); %mess-o-constants 2
    
    %********************************************************
    %****  Compute apply warm layer  correction *************
    %********************************************************
    
    intime=yd-fix(yd);
    loc=(lonx+7.5)/15;
    chktime=loc+intime*24;
    if chktime>24
        chktime=chktime-24;
    end
    newtime=(chktime-24*fix(chktime/24))*3600;
    if icount>1                                  %not first time thru
        if newtime<=21600 | jump==0
            jump=0;
            if newtime < jtime		%re-zero at midnight
                jamset=0;
                fxp=.5;
                tk_pwp=max_pwp;
                tau_ac=0;
                qcol_ac=0;
                dt_wrm=0;
            else
                %************************************
                %****   set warm layer constants  ***
                %************************************
                dtime=newtime-jtime;             %delta time for integrals
                qr_out=Rnl+hs_old+hl_old+RF_old; %total cooling at surface
                q_pwp=fxp*Rns-qr_out;            %tot heat abs in warm layer
                  qqrx(ibg)=hs_old;
                if q_pwp>=50 | jamset==1         %Check for threshold
                    jamset=1;			         %indicates threshold crossed
                    tau_ac=tau_ac+max(.002,tau_old)*dtime;	%momentum integral
                    if qcol_ac+q_pwp*dtime>0	            %check threshold for warm layer existence
                        %******************************************
                        % Compute the absorption profile
                        %******************************************
                        for i=1:5                           %loop 5 times for fxp
                            fxp=1-(0.28*0.014*(1-exp(-tk_pwp/0.014))+0.27*0.357*(1-exp(-tk_pwp/0.357))+0.45*12.82*(1-exp(-tk_pwp/12.82)))/tk_pwp;
                            qjoule=(fxp*Rns-qr_out)*dtime;
                            if qcol_ac+qjoule>0         %Compute warm-layer depth
                                tk_pwp=min(max_pwp,ctd1*tau_ac/sqrt(qcol_ac+qjoule));
                            end;
                        end;
                    else             %warm layer wiped out
                        fxp=0.75;
                        tk_pwp=max_pwp;
                        qjoule=(fxp*Rns-qr_out)*dtime;
                    end;
                    qcol_ac=qcol_ac+qjoule; %heat integral
                    %*******  compute dt_warm  ******
                    if qcol_ac>0
                        dt_wrm=ctd2*(qcol_ac)^1.5/tau_ac;
                    else
                        dt_wrm=0;
                    end;
                end;%                    end threshold check
            end;%                        end midnight reset
            if tk_pwp<ts_depth           %Compute warm layer correction
                dsea=dt_wrm;
            else
                dsea=dt_wrm*ts_depth/tk_pwp;
            end;
        end;%                                    end 6am start first time thru
    end;%                                        end first time thru check
    jtime=newtime;
    %************* output from routine  *****************************
%Output:   A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq  Cd  Ch  Ce  L zet dterx dqerx tkt Urf Trf Qrf RHrf UrfN Rnl  Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis];
%              1   2   3   4   5    6     7    8   9  10  11   12  13  14  15 16  17   18   19    20  21  22  23  24   25   26  27  28  29  30  31    32     33     34   35 36  37   38  39  40  41  42  43 
    ts=tsea+dsea;
    Bx=coare36vn_zrf(u,zu(ibg),t,zt(ibg),rh,zq(ibg),P,ts,Rs,Rl,latx,zi,rain,Ss(ibg),cp(ibg),sigH(ibg),zrf_u(ibg),zrf_t(ibg),zrf_q(ibg));
    tau_old=Bx(2);            %hold stress
    hs_old=Bx(3);             %hold shf
    hl_old=Bx(4);             %hold lhf - use Webb corrected value
    dter=Bx(17);
    RF_old=Bx(34);            %hold rain flux
    
    B(ibg,1)=dt_wrm;   % warming across entire warm layer deg.C
    B(ibg,2)=tk_pwp;   % warm layer thickness m
    B(ibg,3)=dsea;     % heating at selected depth
    
    icount=icount+1;
    
end; %  data line loop

%**************************************************
% Recompute fluxes with warm layer
%**************************************************
clear Bx
Tsea=Tsea+B(:,3)';
Bx=coare36vn_zrf(Ur,zu,Tair,zt,RH,zq,Pair,Tsea,Solar,IR,Lat,zi,Rainrate,Ss,cp,sigH,zrf_u,zrf_t,zrf_q);

B=[Bx B];    %Add the warm layer variables    

%************* output from routine  *****************************
%Output:   A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq  Cd  Ch  Ce  L zet dterx dqerx tkt Urf Trf Qrf RHrf UrfN Rnl  Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis];
%              1   2   3   4   5    6     7    8   9  10  11   12  13  14  15 16  17   18   19    20  21  22  23  24   25   26  27  28  29  30  31    32     33     34   35 36  37   38  39  40  41  42  43 
%       dt_warm tk_pwp dsea    
%          44     45    46