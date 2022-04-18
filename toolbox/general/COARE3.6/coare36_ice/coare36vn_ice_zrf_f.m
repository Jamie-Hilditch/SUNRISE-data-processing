function A=coare36vn_zrf(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q)
%
% Vectorized version of COARE 3.0 code (Fairall et al, 2003) with 
% modification based on the CLIMODE, MBL and CBLAST experiments 
% (Edson et al., 2012). The cool skin option is retained but warm layer 
% and surface wave options removed. 
%
% This version include parameterizations using wave height and wave
% slope using cp and sigH.  If these are set to NaN, then the wind
% speed dependent formulation is used.  The parameterizations are based
%on fits to the Banner-Norison wave model and the Fairall-Edson flux
%database.  It also allows salinity as a input.  Open ocean Ss=35; Great
%Lakes Ss=0;
%
%********************************************************************
% An important component of this code is whether the inputed ts 
% represents the skin temperature of a near surface temperature.  
% How this variable is treated is determined by the jcool parameter:
% set jcool=1 if Ts is bulk ocean temperature (default),
%     jcool=0 if Ts is true ocean skin temperature. 
%********************************************************************

jcool=1;

% The code assumes u,t,rh,ts are vectors; 
% sensor heights zu,zt,zl, latitude lat, and PBL height zi are constants;
% air pressure P and radiation Rs,Rl may be vectors or constants. 
% Default values are assigned for P,Rs,Rl,lat,and zi if these data are not 
% available.  Input NaNs to indicate no data. Defaults should be set to 
% representative regional values if possible.
%
% Input:  
%
%     u = relative wind speed (m/s) at height zu(m)
%     t = bulk air temperature (degC) at height zt(m)
%    rh = relative humidity (%) at height zq(m)
%     P = surface air pressure (mb) (default = 1015)
%    ts = water temperature (degC) see jcool below  if ts<Tfreezing assumes ts=tice
%    Rs = downward shortwave radiation (W/m^2) (default = 150) 
%    Rl = downward longwave radiation (W/m^2) (default = 370)
%   lat = latitude (default = +45 N)
%    zi = PBL height (m) (default = 600m)
%  rain = rain rate (mm/hr)
%    Ss = sea surface salinity (PSU)
%    cp = phase speed of dominant waves (m/s)  
%  sigH =  significant wave height (m)
%  zu, zt, zq heights of the observations (m)
%  zrf_u, zrf_t, zrf_q  reference height for profile  
%
% The user controls the output.  This is currently set as:
%
% Output:  A=[usr tau hsb hlb hbb hsbb tsr qsr zot zoq Cd Ch Ce  L zet dter dqer tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10];
%              1   2   3   4   5   6    7   8   9   10  1  2 3   4  5   6    7    8   9   20  1   2  3   4   5   6   7   8   9  30
%  where
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
%Qs     = surface specific humidity (g/kg)
%Evap   = evaporation rate (mm/h)
%T10    = air temperature at 10m
%Q10    = air specific humidity at 10m
%RH10   = air relative humidity at 10m
%uq     = gustiness velocity (m/s)
%Whf    = whitecap fraction
%Edis   = energy dissipated by wave breaking (W/m^2)

% Notes: 1) u is the relative wind speed, i.e., the magnitude of the
%           difference between the wind (at zu) and ocean surface current 
%           vectors.
%        2) Set jcool=0 in code if ts is true surface skin temperature,
%           otherwise ts is assumed the bulk temperature and jcool=1.
%        3) Set P=NaN to assign default value if no air pressure data 
%           available. 
%        4) Set Rs=NaN, Rl=NaN if no radiation data available.  This assigns 
%           default values to Rs, Rl so that cool skin option can be applied. 
%        5) Set lat=NaN and/or zi=NaN to assign default values if latitude
%           and/or PBL height not given. 
%        6) The code to compute the heat flux caused by precipitation is 
%           included if rain data is available (default is no rain).
%        7) Code updates the cool-skin temperature depression dter and thickness
%           tkt during iteration loop for consistency.
%        8) Number of iterations set to nits = 10.

% Reference:
%
%  Fairall, C.W., E.F. Bradley, J.E. Hare, A.A. Grachev, and J.B. Edson (2003),
%  Bulk parameterization of air sea fluxes: updates and verification for the 
%  COARE algorithm, J. Climate, 16, 571-590.

% Code history:
% 
% 1. 12/14/05 - created based on scalar version coare26sn.m with input
%    on vectorization from C. Moffat.  
% 2. 12/21/05 - sign error in psiu_26 corrected, and code added to use variable
%    values from the first pass through the iteration loop for the stable case
%    with very thin M-O length relative to zu (zetu>50) (as is done in the 
%    scalar coare26sn and COARE3 codes).
% 3. 7/26/11 - S = dt was corrected to read S = ut.
% 4. 7/28/11 - modification to roughness length parameterizations based 
%    on the CLIMODE, MBL, Gasex and CBLAST experiments are incorporated
% 5. New wave parameterization added 9/20/2017  based on fits to wave model
%-----------------------------------------------------------------------
display('flockewe');
% convert input to column vectors
u=u(:);t=t(:);rh=rh(:);P=P(:);ts=ts(:);
Rs=Rs(:);Rl=Rl(:);lat=lat(:);zi=zi(:);
zu=zu(:);zt=zt(:);zq=zq(:);
Ss=Ss(:);cp=cp(:);sigH=sigH(:);
rain=rain(:);
N=length(u);

% set local variables to default values if input is NaN
if isnan(P); P=1013*ones(N,1); end;      % pressure
if isnan(Rs); Rs=150*ones(N,1); end;     % incident shortwave radiation
if isnan(Rl); Rl=370*ones(N,1); end;     % incident longwave radiation
if isnan(lat); lat=45; end;              % latitude
if isnan(zi); zi=600; end;               % PBL height


% input variable u is assumed relative wind speed (magnitude of difference
% between wind and surface current vectors). to follow orginal Fairall code, set
% surface current speed us=0. if us data are available, construct u prior to
% using this code.
us = 0*u;

% convert rh to specific humidity
Tf=-0.0575*Ss+1.71052E-3*Ss.^1.5-2.154996E-4*Ss.*Ss;%freezing point of seawater
Qs = qsat26sea(ts,P,Ss,Tf)./1000;                   % surface specific humidity (g/kg)
Qsw=qsat26sea(Tf,P,Ss,Tf)./1000;                    %Qs over water at freezing pt 
ii=find(ts>Tf);Qsw(ii)=Qs(ii);
[Q,Pv]  = qsat26air(t,P,rh);                        % specific humidity of air (g/kg).  Assumes rh relative to ice T<0
Q=Q./1000;
tsw=ts;ii=find(ts<Tf);tsw(ii)=Tf(ii);               %tsw is water temperature
% 
zoi=5E-4*ones(N,1);                                 %assumed roughness of ice   
%***********  set constants **********************************************
zref=10;
Beta = 1.2;
von  = 0.4;
fdg  = 1.00; % Turbulent Prandtl number
tdk  = 273.16;
grav = grv(lat);

%***********  air constants **********************************************
Rgas = 287.1;
Le   = (2.501-.00237*ts)*1e6;
cpa  = 1004.67;
cpv  = cpa*(1+0.84*Q);
rhoa = P*100./(Rgas*(t+tdk).*(1+0.61*Q));
rhodry = (P-Pv)*100./(Rgas*(t+tdk));
visa = 1.326e-5*(1+6.542e-3.*t+8.301e-6*t.^2-4.84e-9*t.^3);

%***********  cool skin constants  ***************************************
Al35   = 2.1e-5*(tsw+3.2).^0.79;
Al0   =(2.2*real((tsw-1).^0.82)-4)*1e-5;
Al=Al0+(Al35-Al0).*Ss/35;
bets=7.5e-4;%assumes beta a constant
be   = bets*Ss;%be is beta*Salinity
cpw  = 4000;
rhow = 1022;
visw = 1e-6;
tcw  = 0.6;
bigc = 16*grav*cpw*(rhow*visw)^3./(tcw.^2*rhoa.^2);
wetc = 0.622*Le.*Qs./(Rgas*(tsw+tdk).^2);

%***********  net radiation fluxes ***************************************
Rns = 0.945.*Rs; % albedo correction
% IRup = eps*sigma*T^4 + (1-eps)*IR
% Rnl = IRup - IR
% Rnl = eps*sigma*T^4 - eps*IR  as below

Rnl = 0.97*(5.67e-8*(ts-0.3*jcool+tdk).^4-Rl); % initial value

% IRup = Rnl + IR

%****************  begin bulk loop ********************************************

%***********  first guess ************************************************
du = u-us;
dt = ts-t-.0098.*zt;
dtw=dt;ii=find(ts<Tf);dtw(ii)=Tf(ii)-t(ii)-.0098.*zt(ii);%sea-air temp diff over water
dq = Qs-Q;
dqw=dq;ii=find(ts<Tf);dqw(ii)=Qsw(ii)-Q(ii);%sea-air q diff over water
ta = t+tdk;
ug = 0.5;
dter  = 0.3;
ut    = sqrt(du.^2+ug.^2);
u10   = ut.*log(10/1e-4)./log(zu/1e-4);
usr   = 0.035*u10;
zo10  = 0.011*usr.^2./grav + 0.11*visa./usr;
Cd10  = (von./log(10./zo10)).^2;
Ch10  = 0.00115;
Ct10  = Ch10./sqrt(Cd10);
zot10 = 10./exp(von./Ct10);
Cd    = (von./log(zu./zo10)).^2;
Ct    = von./log(zt./zot10);
CC    = von*Ct./Cd;
Ribcu = -zu./zi./.004/Beta^3;
Ribu  = -grav.*zu./ta.*((dt-dter*jcool)+.61*ta.*dq)./ut.^2;
zetu = CC.*Ribu.*(1+27/9*Ribu./CC);
k50=find(zetu>50); % stable with very thin M-O length relative to zu
k=find(Ribu<0); 
if length(Ribcu)==1
    zetu(k)=CC(k).*Ribu(k)./(1+Ribu(k)./Ribcu); clear k;
else
    zetu(k)=CC(k).*Ribu(k)./(1+Ribu(k)./Ribcu(k)); clear k;
end
L10 = zu./zetu;
gf=ut./du;
usrw = ut.*von./(log(zu./zo10)-psiu_40(zu./L10));
tsrw = -(dtw-dter*jcool).*von*fdg./(log(zt./zot10)-psit_26(zt./L10));
qsrw = -(dqw-wetc.*dter*jcool)*von*fdg./(log(zq./zot10)-psit_26(zq./L10));
usri = ut.*von./(log(zu./zo10)-psiu_40(zu./L10));
tsri = -(dt).*von*fdg./(log(zt./zot10)-psit_26(zt./L10));
qsri = -(dq)*von*fdg./(log(zq./zot10)-psit_26(zq./L10));

tkt = 0.001*ones(N,1);

%**********************************************************
%  The following gives the new formulation for the
%  Charnock variable
%**********************************************************
%%%%%%%%%%%%%   COARE 3.5 wind speed dependent charnock
charnC = 0.011*ones(N,1);
umax=19;
a1=0.0017;
a2=-0.0050;
charnC=a1*u10+a2;
k=find(u10>umax);
charnC(k)=a1*umax+a2;


%%%%%%%%%   if wave age is given but not wave height, use parameterized
%%%%%%%%%   wave height
hsig=(0.02*(cp./u10).^1.1-0.0025).*u10.^2;
hsig=max(hsig,.25);
ii=find(~isnan(cp) & isnan(sigH));
sigH(ii)=hsig(ii);

Ad=0.15;  %Sea-state/wave-age dependent coefficients from wave model
Bd=2.2;
zoS=sigH.*Ad.*(usr./cp).^Bd;
charnS=zoS.*grav./usr./usr;

nits=10; % number of iterations
charn=charnC;
ii=find(~isnan(cp));charn(ii)=charnS(ii);
%**************  bulk loop **************************************************
% usrw=usr;usri=usr;
% tsrw=tsr;tsri=tsr;
% qsrw=qsr;qsri=qsr;
utw=ut;uti=ut;
for i=1:nits
    
    zetw=von.*grav.*zu./ta.*(tsrw +.61*ta.*qsrw)./(usrw.^2);
    Lw=zu./zetw;
    zeti=von.*grav.*zu./ta.*(tsri +.61*ta.*qsri)./(usri.^2);
    Li=zu./zeti;
    zow=charn.*usrw.^2./grav+0.11*visa./usrw; % surface roughness
    zoI=zoi+0.11*visa./usri;
    rrw=zow.*usrw./visa;
    rri=zoI.*usri./visa;
    zoqw=min(1.6e-4,5.8e-5./rrw.^.72);       % These thermal roughness lengths give Stanton and
    zotw=zoqw;     
         
    %%%%%%%%%%%%%  %Andreas 1987 for snow/ice
    
   		rti=rri*exp(1.250);
     	rqi=rri*exp(1.610);
    ik=find(rri>.135 & rri<=2.5);
        rti(ik)=rri(ik).*exp(0.149-.55*log(rri(ik)));
     	rqi(ik)=rri(ik).*exp(0.351-0.628*log(rri(ik)));
     ik=find(rri>2.5 & rri<=1000);
     	rti(ik)=rri(ik).*exp(0.317-0.565*log(rri(ik))-0.183*log(rri(ik)).*log(rri(ik)));
      	rqi(ik)=rri(ik).*exp(0.396-0.512*log(rri(ik))-0.180*log(rri(ik)).*log(rri(ik)));

      zoti=rti.*visa./usri; zoqi=rqi.*visa./usri;
     
    cdhfw=von./(log(zu./zow)-psiu_26(zu./Lw));
    cqhfw=von.*fdg./(log(zq./zoqw)-psit_26(zq./Lw));
    cthfw=von.*fdg./(log(zt./zotw)-psit_26(zt./Lw));
    usrw=utw.*cdhfw;
    qsrw=-(dqw-wetc.*dter*jcool).*cqhfw;
    tsr=-(dtw-dter*jcool).*cthfw;
    tvsrw=tsrw+0.61*ta.*qsrw;
    tssrw=tsrw+0.51*ta.*qsrw;
    Bfw=-grav./ta.*usrw.*tvsrw;
    ugw=0.2*ones(N,1);
    k=find(Bfw>0); 
    if length(zi)==1;
        ugw(k)=Beta*(Bfw(k).*zi).^.333; clear k;
    else
        ugw(k)=Beta*(Bfw(k).*zi(k)).^.333; clear k;
    end
    utw=sqrt(du.^2+ugw.^2);
    gfw=utw./du;
    hsbw=-rhoa*cpa.*usrw.*tsrw;
    hlbw=-rhoa.*Le.*usrw.*qsrw;
    
    cdhfi=von./(log(zu./zoI)-psiu_26(zu./Li));
    cqhfi=von.*fdg./(log(zq./zoqi)-psit_26(zq./Li));
    cthfi=von.*fdg./(log(zt./zoti)-psit_26(zt./Li));
    usri=uti.*cdhfi;
    qsri=-(dq).*cqhfi;
    tsri=-(dt).*cthfi;
    tvsri=tsri+0.61*ta.*qsri;
    tssri=tsri+0.51*ta.*qsri;
    Bfi=-grav./ta.*usri.*tvsri;
    ugi=0.2*ones(N,1);
    k=find(Bfi>0); 
    if length(zi)==1;
        ugi(k)=Beta*(Bfi(k).*zi).^.333; clear k;
    else
        ugi(k)=Beta*(Bfi(k).*zi(k)).^.333; clear k;
    end
    uti=sqrt(du.^2+ugi.^2);
    gfi=uti./du;
    hsbi=-rhoa*cpa.*usri.*tsri;
    hlbi=-rhoa.*Le.*usri.*qsri;
    
    qout=Rnl+hsbw+hlbw;
    dels=Rns.*(0.065+11*tkt-6.6e-5./tkt.*(1-exp(-tkt/8.0e-4)));
    qcol=qout-dels;
    alq=Al.*qcol+be.*hlbw.*cpw./Le;
    xlamx=6.0*ones(N,1);
    tkt=min(0.01, xlamx.*visw./(sqrt(rhoa./rhow).*usrw));
    k=find(alq>0); xlamx(k)=6./(1+(bigc(k).*alq(k)./usrw(k).^4).^0.75).^0.333;
    tkt(k)=xlamx(k).*visw./(sqrt(rhoa(k)./rhow).*usrw(k)); clear k;
    dter=qcol.*tkt./tcw;
    dqer=wetc.*dter;
    Rnl=0.97*(5.67e-8*(ts-dter.*jcool+tdk).^4-Rl); % update dter
    if i==1; % save first iteration solution for case of zetu>50;
        usrw50=usrw(k50);tsrw50=tsrw(k50);qsrw50=qsrw(k50);Lw50=Lw(k50);
        zetw50=zetw(k50);dter50=dter(k50);dqer50=dqer(k50);tkt50=tkt(k50);
    end
    u10Nw = usrw./von./gfw.*log(10./zow);
    charnC=a1*u10Nw+a2;
    k=find(u10Nw>umax);
    charnC(k)=a1*umax+a2;
    charn=charnC;
    zoS=sigH.*Ad.*(usr./cp).^Bd;%-0.11*visa./usr;
    charnS=zoS.*grav./usr./usr;
    ii=find(~isnan(cp));charn(ii)=charnS(ii);
end

% insert first iteration solution for case with zetu>50
usrw(k50)=usrw50;tsrw(k50)=tsrw50;qsrw(k50)=qsrw50;Lw(k50)=Lw50;
zetw(k50)=zetw50;dter(k50)=dter50;dqer(k50)=dqer50;tkt(k50)=tkt50;

%****************  compute WATER fluxes  ********************************************
tauw=rhoa.*usrw.*usrw./gfw;      % wind stress
hsbw=-rhoa.*cpa.*usrw.*tsrw;     % sensible heat flux
hlbw=-rhoa.*Le.*usrw.*qsrw;      % latent heat flux
hbbw=-rhoa.*cpa.*usrw.*tvsrw;    % buoyancy flux
hsbbw=-rhoa.*cpa.*usrw.*tssrw;   % sonic heat flux
wbarw=1.61*hlbw./Le./(1+1.61*Q)./rhoa+hsbw./rhoa./cpa./ta;
hlwebbw=hlbw+rhoa.*wbarw.*Q.*Le;
Evapw=1000*hlwebbw./Le./1000*3600;   %mm/hour
%****************  compute cie fluxes  ********************************************
taui=rhoa.*usri.*usri./gfi;      % wind stress
hsbi=-rhoa.*cpa.*usri.*tsri;     % sensible heat flux
hlbi=-rhoa.*Le.*usri.*qsri;      % latent heat flux
hbbi=-rhoa.*cpa.*usri.*tvsri;    % buoyancy flux
hsbbi=-rhoa.*cpa.*usri.*tssri;   % sonic heat flux
wbari=1.61*hlbi./Le./(1+1.61*Q)./rhoa+hsbi./rhoa./cpa./ta;
hlwebbi=hlbi+rhoa.*wbari.*Q.*Le;
Evapi=1000*hlwebbi./Le./1000*3600;   %mm/hour

%*****  compute transfer coeffs relative to ut @ meas. ht  ********************
Cdw=tauw./rhoa./utw./max(.1,du);
Chw=-usrw.*tsrw./utw./(dtw-dter*jcool);
Cew=-usrw.*qsrw./(dqw-dqer*jcool)./ut;
Cdi=taui./rhoa./utw./max(.1,du);
Chi=-usri.*tsri./uti./(dt);
Cei=-usri.*qsri./(dq)./uti;

%***  compute 10-m neutral coeff relative to ut (output if needed) ************
Cdn_10w=1000*von.^2./log(10./zow).^2;
Chn_10w=1000*von.^2.*fdg./log(10./zow)./log(10./zotw);
Cen_10w=1000*von.^2.*fdg./log(10./zow)./log(10./zoqw);

Cdn_10i=1000*von.^2./log(10./zoI).^2;
Chn_10i=1000*von.^2.*fdg./log(10./zoI)./log(10./zoti);
Cen_10i=1000*von.^2.*fdg./log(10./zoI)./log(10./zoqi);


%***  compute 10-m neutral coeff relative to ut (output if needed) ************
%  Find the stability functions
%*********************************
%zrf_u=10;             %User defined reference heights
%zrf_t=10;
%zrf_q=10;
psiw=psiu_26(zu./Lw);
psi10w=psiu_26(10./Lw);
psirfw=psiu_26(zrf_u./Lw);
psiTw=psit_26(zt./Lw);
psi10Tw=psit_26(10./Lw);
psirfTw=psit_26(zrf_t./Lw);
psirfQw=psit_26(zrf_q./Lw);
gfw=utw./du;

%*********************************************************
%  Determine the wind speeds relative to ocean surface
%  Note that usr is the friction velocity that includes 
%  gustiness usr = sqrt(Cd) S, which is equation (18) in
%  Fairall et al. (1996)
%*********************************************************
Sw = utw;
U = du;
S10w = Sw + usrw./von.*(log(10./zu)-psi10w+psiw);
U10w = S10w./gfw;
% or U10 = U + usr./von./gf.*(log(10/zu)-psi10+psi);
Urfw = U + usrw./von./gfw.*(log(zrf_u./zu)-psirfw+psiw);
UNw = U + psiw.*usrw/von./gfw;
U10Nw = U10w + psi10w.*usrw/von./gfw;
UrfNw = Urfw + psirfw.*usrw/von./gfw;

UN2w = usrw/von./gfw.*log(zu./zow);
U10N2w = usrw./von./gfw.*log(10./zow);
UrfN2w  = usrw./von./gfw.*log(zrf_u./zow);

%******** rain heat flux (save to use if desired) *****************************
if isnan(rain(1));
    RF=zeros(size(usrw));
else
    dwat=2.11e-5*((t+tdk)./tdk).^1.94; %! water vapour diffusivity
    dtmp=(1. + 3.309e-3*t - 1.44e-6.*t.*t).*0.02411./(rhoa.*cpa); %! heat diffusivity
    dqs_dt=Q.*Le./(Rgas.*(t+tdk).^2); %! Clausius-Clapeyron
    alfac= 1./(1+0.622*(dqs_dt.*Le.*dwat)./(cpa.*dtmp)); %! wet bulb factor
    RF= rain.*alfac.*cpw.*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool).*Le./cpa)./3600;
end

lapse=grav/cpa;
SST=ts-dter*jcool;

T = t;
%[size(-psi10T+psiT + lapse*(zt-10))]
T10w = T + tsrw./von.*(log(10./zt)-psi10Tw+psiTw) + lapse.*(zt-10);
Trfw = T + tsrw./von.*(log(zrf_t./zt)-psirfTw+psiTw) + lapse.*(zt-zrf_t);
TNw = T + psiTw.*tsrw/von;
T10Nw = T10w + psi10Tw.*tsrw/von;
TrfNw = Trfw + psirfTw.*tsrw/von;

TN2w = SST + tsrw/von.*log(zt./zotw)-lapse.*zt;
T10N2w = SST + tsrw/von.*log(10./zotw)-lapse.*10;
TrfN2w = SST + tsrw/von.*log(zrf_t./zotw)-lapse.*zrf_t;

dqer=wetc.*dter*jcool;
SSQ=Qs-dqer;
SSQ=SSQ*1000;
Q=Q*1000;
qsrw=qsrw*1000;
Q10w = Q + qsrw./von.*(log(10./zq)-psi10Tw+psiTw);
Qrfw = Q + qsrw./von.*(log(zrf_q./zq)-psirfQw+psiTw);
QNw = Q + psiTw.*qsrw/von./sqrt(gfw);
Q10Nw = Q10w + psi10Tw.*qsrw/von;
QrfNw = Qrfw + psirfQw.*qsrw/von;

QN2w = SSQ + qsrw/von.*log(zq./zoqw);
Q10N2w = SSQ + qsrw/von.*log(10./zoqw);
QrfN2w = SSQ + qsrw/von.*log(zrf_q./zoqw);
RHrfw=RHcalc(Trfw,P,Qrfw/1000,Tf);
RH10w=RHcalc(T10w,P,Q10w/1000,Tf);

%%%%%%%%%%%%%%%%%  Other wave breaking statistics fromn Banner-Morison wave
%%%%%%%%%%%%%%%%%  model
Whf=1.25e-3*U10w.^1.1./sqrt(cp./U10w);%  whitecap fraction
Edis=0.095*rhoa.*U10w.*usrw.^2;             %  energy dissipation rate from breaking waves W/m^2


%****************  output  ****************************************************
usr=usrw;tau=tauw;hsb=hsbw;hlb=hlbw;hbb=hbbw;hsbb=hsbbw;hlwebb=hlwebbw;tsr=tsrw;qsr=qsrw;zo=zow;zot=zotw;zoq=zoqw;Cd=Cdw;Ch=Chw;Ce=Cew;L=Lw;zet=zetw;Urf=Urfw;Trf=Trfw;Qrf=Qrfw;RHrf=RHrfw;UrfN=UrfNw;
UN=UNw;U10=U10w;U10N=U10Nw;Cdn_10=Cdn_10w;Chn_10=Chn_10w;Cen_10=Cen_10w;Evap=Evapw;T10=T10w;Q10=Q10w;RH10=RH10w;ug=ugw;
A=[usr tau hsb hlb hbb hsbb hlwebb tsr qsr zo  zot zoq Cd Ch Ce  L zet dter dqer tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis usri taui hsbi hlbi hbbi hsbbi zeti Cdn_10i Chn_10i Cen_10i Cdi Chi Cei];
%   1   2   3   4   5   6    7      8   9  10  11  12  13 14 15 16  17   18   19  20  21  22  23   24   25 26  27  28  29  30    31     32     33   34 35  36  37  38   39  40  41  42  43   44    45  46    47  48   49    50     51
end
%------------------------------------------------------------------------------
function psi=psit_26(zet)
% computes temperature structure function
dzet=min(50,0.35*zet); % stable
psi=-((1+0.6667*zet).^1.5+0.6667*(zet-14.28).*exp(-dzet)+8.525);
k=find(zet<0); % unstable
x=(1-15*zet(k)).^0.5;
psik=2*log((1+x)./2);
x=(1-34.15*zet(k)).^0.3333;
psic=1.5*log((1+x+x.^2)./3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
f=zet(k).^2./(1+zet(k).^2);
psi(k)=(1-f).*psik+f.*psic;
end
%------------------------------------------------------------------------------
function psi=psiu_26(zet)
% computes velocity structure function
dzet=min(50,0.35*zet); % stable
a=0.7;
b=3/4;
c=5;
d=0.35;
psi=-(a*zet+b*(zet-c/d).*exp(-dzet)+b*c/d);
k=find(zet<0); % unstable
x=(1-15*zet(k)).^0.25;
psik=2*log((1+x)/2)+log((1+x.*x)/2)-2*atan(x)+2*atan(1);
x=(1-10.15*zet(k)).^0.3333;
psic=1.5*log((1+x+x.^2)/3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
f=zet(k).^2./(1+zet(k).^2);
psi(k)=(1-f).*psik+f.*psic;
end
%------------------------------------------------------------------------------
function psi=psiu_40(zet)
% computes velocity structure function
dzet=min(50,0.35*zet); % stable
a=1;
b=3/4;
c=5;
d=0.35;
psi=-(a*zet+b*(zet-c/d).*exp(-dzet)+b*c/d);
k=find(zet<0); % unstable
x=(1-18*zet(k)).^0.25;
psik=2*log((1+x)/2)+log((1+x.*x)/2)-2*atan(x)+2*atan(1);
x=(1-10*zet(k)).^0.3333;
psic=1.5*log((1+x+x.^2)/3)-sqrt(3)*atan((1+2*x)./sqrt(3))+4*atan(1)./sqrt(3);
f=zet(k).^2./(1+zet(k).^2);
psi(k)=(1-f).*psik+f.*psic;
end
%------------------------------------------------------------------------------
function exx=bucksat(T,P,Tf)
% computes saturation vapor pressure [mb]
% given T [degC] and P [mb] Tf is freezing pt 
exx=6.1121.*exp(17.502.*T./(T+240.97)).*(1.0007+3.46e-6.*P);
ii=find(T<Tf);
exx(ii)=(1.0003+4.18e-6*P(ii)).*6.1115.*exp(22.452.*T(ii)./(T(ii)+272.55));%vapor pressure ice
end
%------------------------------------------------------------------------------
function qs=qsat26sea(T,P,Ss,Tf)
% computes surface saturation specific humidity [g/kg]
% given T [degC] and P [mb]
ex=bucksat(T,P,Tf);
fs=1-0.02*Ss/35;% reduction sea surface vapor pressure by salinity
es=fs.*ex; 
qs=622*es./(P-0.378*es);
end
%------------------------------------------------------------------------------
function [q,em]=qsat26air(T,P,rh)
% computes saturation specific humidity [g/kg]
% given T [degC] and P [mb]
Tf=0;%assumes relative humidity for pure water
es=bucksat(T,P,Tf);
em=0.01*rh.*es;
q=622*em./(P-0.378*em);
end
%------------------------------------------------------------------------------
function g=grv(lat)
% computes g [m/sec^2] given lat in deg
gamma=9.7803267715;
c1=0.0052790414;
c2=0.0000232718;
c3=0.0000001262;
c4=0.0000000007;
phi=lat*pi/180;
x=sin(phi);
g=gamma*(1+c1*x.^2+c2*x.^4+c3*x.^6+c4*x.^8);
end
%------------------------------------------------------------------------------
function RHrf=RHcalc(T,P,Q,Tf)
% computes relative humidity given T,P, & Q
es=6.1121.*exp(17.502.*T./(T+240.97)).*(1.0007+3.46e-6.*P);
ii=find(T<Tf);%ice case
es(ii)=6.1115.*exp(22.452.*T(ii)./(T(ii)+272.55)).*(1.0003+4.18e-6*P(ii));
em=Q.*P./(0.378.*Q+0.622);
RHrf=100*em./es;
end