disp('metz_red5');

load('C:\data\cwf\Dropbox\DATA\HIWINGS_2013\Knorr\flux\flux_HIWINGS');;%read file with hour average data; set your local directory 

flux2=flux3;
jdy=flux2(:,1);%day-of-year at beginning of time average
ushp=flux2(:,2);%doppler log, SCS (m/s) ?not available, replaced by SOG
U=flux2(:,3);%true wind PSD sonic relative to water (m/s)
dir=flux2(:,4);%true wind direction from relative to water, PSD sonic (deg)
urel=flux2(:,5);%relative wind speed, PSD (m/s)
reldir=flux2(:,6);%relative wind dir (from),clockwise rel ship's bow, PSD sonic (deg)
head=flux2(:,7);%ship heading, deg clockwise rel north, PSD GPS
tsnk=flux2(:,8);%sea snake temperature, PSD (C)
ta=flux2(:,9);%air temperature, PSD (C)%
qse=flux2(:,10);%sea surface specific humidity from snake (g/kg)
qa=flux2(:,11);%air specific humidity, PSD (g/kg)
hsc=flux2(:,12);%sensible heat flux, covariance, PSD sonic anemometer(W/m^2)
hsib=flux2(:,13);%sensible heat flux, ID, PSD sonic anemometer(W/m^2)
hsb=flux2(:,14);%bulk sensible heat flux, (W/m^2)
hlc=flux2(:,15);%latent heat flux, covariance, (W/m^2)
hlib=flux2(:,16);%latent heat flux, ID, (W/m^2)
hlb=flux2(:,17);%bulk latent heat flux, W/m^2 (includes Webb et al. correction)
taucx=flux2(:,18);%covariance streamwise stress, PSD sonic anemometer (N/m^2)
taucy=flux2(:,19);%covariance cross-stream stress, PSD sonic anemometer (N/m^2)
tauib=flux2(:,20);%ID streamwise stress, PSD sonic anemometer (N/m^2)
taub=flux2(:,21);%bulk wind stress along mean wind, (N/m^2)
rs=flux2(:,22);%downward solar flux, PSD units (W/m^2)
rl=flux2(:,23);%downward IR flux, PSD units (W/m^2)
org=flux2(:,24);%rainrate, PSD STI optical rain gauge (mm/hr)
J=flux2(:,25);%ship plume contamination index
tiltx=flux2(:,26);%flow tilt at PSD sonic anemometer
Jm=flux2(:,27);%ship maneuver index
ct=flux2(:,28);%ct^2 (K^2/m^.667)
cq=flux2(:,29);%cq^2 ((g/kg)^2/m^.667)
cu=flux2(:,30);%cu^2 ((m/s)^2/m^.667)
cw=flux2(:,31);%cw^2 ((m/s)^2/m^.667)
hrain=flux2(:,32);%rain heat flux, Gosnell et al 1995, JGR, 18437-18442 (W/m^2)
hlwebb=flux2(:,33);%correction to measured latent heat flux, Webb et al.
lat=flux2(:,34);%latitude, deg (PSD pcode)
lon=flux2(:,35);%longitude, deg (PSD pcode)
zu_etl=flux2(:,36);%height of mean wind sensor, 17.75 m
zt_etl=flux2(:,37);%height of mean air temperature sensor, 15.6 m
zq_etl=flux2(:,38);%height of mean air humidity sensor, 15.6 m
%***** ships imet and scs data
sog2=flux2(:,39);%speed over ground, SCS gps, (m/s)
U_scs=flux2(:,40); %true wind speed relative to earth m/s) – 19.6m
dir_scs=flux2(:,41);%true wind direction from relative to earth, clockwise rel north (deg)
cog2=flux2(:,42);%%course over ground, SCS gps, (m/s)
tsg=flux2(:,43);%tsg water temperature (C)
ta_im=flux2(:,44);%imet air temperature (C)
qs_tsg=flux2(:,45);%imet bulk water specific humidity (g/kg)
qa_im=flux2(:,46);%imet air specific humidity, (g/kg)
rs_im=flux2(:,47);%imet solar flux, (W/m^2)
rl_im=flux2(:,48);%imet IR flux (W/m^2)
wco2_lic=flux2(:,49);%LICOR 7500 CO2 flux, (micatm m/s)
q_lic=flux2(:,50);%Specific humidity from LICOR (g/kg)
sgq_lic=flux2(:,51);%Standard deviation of specific humidity from LICOR (g/kg)
co2_lic=flux2(:,52); %CO2 concentration from Licor (umol/mol)
sgC_lic=flux2(:,53);%Standard deviation of CO2 concentration from LICOR (microatm)
press=flux2(:,54); %Atmospheric pressure (mb)
Uearth=flux2(:,55); %True wind speed (m/s) relative to earth
direarth=flux2(:,56); % True wind direction (deg) from relative to earth%  ii=find(lat>100);lat(ii)=lat(ii-1);
hent_bulk=flux2(:,57);
hsc_uh=flux2(:,60);%%%%%%%%%%%   uh means the Univ Hawaii system
hsib_uh=flux2(:,61);
taucx_uh=flux2(:,62);
taucy_uh=flux2(:,63);
tauib_uh=flux2(:,64);
rh=relhum5([ta qa press ]);
zu_etl=ones(length(jdy),1)*15.9;%%%%%%%%  sensor height during hiwings
zt_etl=ones(length(jdy),1)*14;
zq_etl=ones(length(jdy),1)*14;

%%%%%%%%%%%  consensus fluxes from PSD, UH, and Ming
taucxm=flux3(:,74);%streamwise stress
tauxym=flux3(:,75);%cross stream stress
tauibm=flux3(:,76);%ID stress
taucc=flux3(:,77);%combines cov and ID
usrcc=flux3(:,78);%friction veloc
hscc=flux3(:,79);%sensible heat flux
%%%%%%%%%%%%%%%%  consensus waves from Reigl, waverider, and WWIII
hsgm=flux3(:,80);%total wave field wave height
cpm=flux3(:,81);%total wave field peak phase speed
hsgms=flux3(:,82);%wind sea wave field
cpms=flux3(:,83);
Wdirm=flux3(:,84);%total wave direction
%%%%%%%%%%%%%%%  Chem variables
dDMS=flux3(:,85);%sea-air delta for dms
k660_dms=flux3(:,86);%transfer coeff at Sc=660 (m/s)
ald=flux3(:,87);%dimensionless solubility
scwd=flux3(:,88);%water side schmidt number
dCO2=flux3(:,89);
k660_co2=flux3(:,90);
alc=flux3(:,91);
scwc=flux3(:,92);
Wf=flux3(:,93);%observed whitecap fraction
Wfv=flux3(:,94);%parameterized whitecap fraction

ffc=(1.58+1.95*log(qa))/10;
ii=find(~isnan(qa));
hlc(ii)=hlc(ii)+(qa(ii)-q_lic(ii))./max(1,qa(ii)).*ffc(ii).*hlc(ii);
hlib(ii)=hlib(ii)+(qa(ii)-q_lic(ii))./max(1,qa(ii)).*ffc(ii).*hlib(ii);
hlib_lic=hlib;
   
%*******************  set constants  ****************
ii=find(isnan(lat));lat(ii)=53;
tdk=273.16;
grav=grv(lat);ii=find(isnan(grav));grav(ii)=9.83;%9.82;
Rgas=287.1;
cpa=1004.67;
rhoair=press.*100./(Rgas*(ta+tdk).*(1+0.61*qa/1000));rhoa=rhoair;
visa=1.326e-5*(1+6.542e-3*ta+8.301e-6*ta.*ta-4.84e-9*ta.*ta.*ta);%kinematic viscocity of air
ii=find(isfinite(tsnk));Le=(2.501-.00237*median(tsnk(ii)))*1e6;

%*******************  set longwave rad flux  ****************
rlclr=(.52+.13/60*abs(lat)+(.082-.03/60*abs(lat)).*sqrt(qa)).*(5.67e-8*(ta+tdk).^4); %clear sky
rlnet=0.97*(rl-5.67e-8*(tsnk+tdk).^4);
rlnetc=0.97*(rlclr-5.67e-8*(tsnk+tdk).^4);

%%%%%%%%%%%%%%%  set variables for flux algorithm inputs
n=length(U);
u=U;           %wind speed 2 to 20 m/s
zu=zu_etl;%ht of wiind speed
t=ta;%air temp (C)
zt=zt_etl;%height of air tem
%rh=80*ones(1,n);%relative humidity
zq=zq_etl;%height of RH
P=press;%pressure (mb)
ts=tsnk;%surface temperature
Rs=rs;%downward solar flux
Rl=rl;%downward IR flux
%lat=20*ones(1,n);%latitude
zi=600*ones(1,n);%mixed layer height
rain=org;%rain rate (mm/hr)
Ss=35*ones(n,1);%salinity (PSU)
visw = SW_Kviscosity2(ts,Ss);
%visw = SW_Kviscosity(tsnk,'C',Ss','ppt');%kinematic viscosity seawater

%%%%%%%%  set reference height for 2nd means 
zrf_u=10;
zrf_t=10;
zrf_q=10;
%%%%%%%%%%%%%%

B=coare35vn(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,cpm*NaN,hsgm*NaN);                              %COARE 35 wind speed
A0=coare36vn_zrf(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,NaN*cpm,NaN*hsgm,zrf_u,zrf_t,zrf_q);    %COARE 36 wind speed
A=coare36vn_zrf(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cpm,hsgm,zrf_u,zrf_t,zrf_q);             %COARE 36 wave  model
Gj=3600*100;
     %COARE36G wind speed
[Adu Gdu] = coareG36vn_dms(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,NaN*cpm,NaN*hsgm,zrf_u,zrf_t,zrf_q,dDMS);    
[Acu Gcu] = coareG36vn_co2(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,NaN*cpm,NaN*hsgm,zrf_u,zrf_t,zrf_q,dCO2);          
k660_co2_Coaru=Gcu(:,14)*Gj;k0660_co2_Coaru=Gcu(:,4).*sqrt(scwc/660)*Gj;kbb660_co2_Coaru=Gcu(:,11).*sqrt(scwc/660)*Gj;
k660_dms_Coaru=Gdu(:,14)*Gj;k0660_dms_Coaru=Gdu(:,4).*sqrt(scwd/660)*Gj;kbb660_dms_Coaru=Gdu(:,11).*sqrt(scwd/660)*Gj;
     %COARE36G wave reynolds number
[Adw Gdw] = coareG36vn_dms(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cpm,hsgm,zrf_u,zrf_t,zrf_q,dDMS);    
[Acw Gcw] = coareG36vn_co2(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cpm,hsgm,zrf_u,zrf_t,zrf_q,dCO2);          
k660_co2_Coarw=Gcw(:,14)*Gj;k0660_co2_Coarw=Gcw(:,4).*sqrt(scwc/660)*Gj;kbb660_co2_Coarw=Gcw(:,11).*sqrt(scwc/660)*Gj;
k660_dms_Coarw=Gdw(:,14)*Gj;k0660_dms_Coarw=Gdw(:,4).*sqrt(scwd/660)*Gj;kbb660_dms_Coarw=Gdw(:,11).*sqrt(scwd/660)*Gj;


taub=B(:,2);hsb=B(:,3);hlb=B(:,4);
wtb=hsb./rhoair./cpa;
wqb=hlb./rhoair./Le*1e3;
wtvb=wtb+.61e-3.*wqb.*(ta+tdk);
wtb_son=wtb+.51e-3.*wqb.*(ta+tdk);
np=length(jdy);
for i=1:np
  	if wtvb(i)>0   
     	wstr(i)=  ( grav(i)*wtvb(i)*600./(ta(i)+tdk))^.333;
  	else
  		wstr(i)=.2/1.2;   
  	end;
end;
   
wgust=1.2*wstr';
gf=sqrt(1+wgust.^2./U.^2);
usb=sqrt(taub./rhoair.*gf);%
tsb=-wtb./usb;
tsbs=-wtb_son./usb;
qsb=-wqb./usb;

hson=hsc+0.51*hlb/Le*cpa.*(ta+tdk);
wbar=1.61*hlb/Le./(1+1.61*qa/1e3)./rhoair+hsb./rhoair/cpa./(ta+tdk);%formulation in hlb already includes webb
hl_webb=rhoair.*wbar.*qa*Le/1000;
ii=find(isnan(hlwebb));
hlwebb(ii)=hl_webb(ii);

zetx=-.4*grav.*zu_etl.*wtvb./usb.^3./(ta+tdk);
Lb=zetx./zu_etl;
zetx=max(zetx,-50);
zetx=min(zetx,50);
u10n=U-usb/.4.*(log(zu_etl/10)-psiu_30(zetx))./gf;
q10n=qa-qsb/.4.*(log(zq_etl/10)-psit_30(zetx.*zq_etl./zu_etl));
th10n=ta-tsb/.4.*(log(zq_etl/10)-psit_30(zetx.*zq_etl./zu_etl))+.0098*zq_etl;

binsize=2;xmin=0;xmax=26;[X, Y, M, SD , N] = binave2(u10n',tauibm',binsize,xmin,xmax);
tauibn=M;u10nbac=X;sigtib=SD;
binsize=2;xmin=0;xmax=26;[X, Y, M, SD,  N] = binave2(u10n',taucxm',binsize,xmin,xmax);
taucxmn=M;u10nbc=X;sigtcx=SD;
binsize=2;xmin=0;xmax=26;[X, Y, M, SD,  N] = binave2(u10n',taub',binsize,xmin,xmax);
taubn=M;u10nbc=X;sigtb=SD;
binsize=2;xmin=0;xmax=26;[X, Y, M, SD , N] = binave2(u10n',rhoair',binsize,xmin,xmax);
rhoairn=M;

binsize=2;xmin=0;xmax=26;[X, Y, M, SD , N] = binave2(u10n',k660_co2_Coaru',binsize,xmin,xmax);
k660_co2_um=M;u10ncu=X;sigcu=SD;
binsize=2;xmin=0;xmax=26;[X, Y, M, SD,  N] = binave2(u10n',k660_co2_Coarw',binsize,xmin,xmax);
k660_co2_wm=M;u10ncw=X;sigcw=SD;
binsize=2;xmin=0;xmax=26;[X, Y, M, SD,  N] = binave2(u10n',k660_co2',binsize,xmin,xmax);
k660_co2_obs=M;u10nco=X;sigco=SD;

binsize=2;xmin=0;xmax=26;[X, Y, M, SD,  N] = binave2(u10n',k660_dms_Coaru',binsize,xmin,xmax);
k660_dms_um=M;u10ndu=X;sigdu=SD;
binsize=2;xmin=0;xmax=26;[X, Y, M, SD , N] = binave2(u10n',k660_dms_Coarw',binsize,xmin,xmax);
k660_dms_wm=M;u10ndw=X;sigdw=SD;
binsize=2;xmin=0;xmax=26;[X, Y, M, SD,  N] = binave2(u10n',k660_dms',binsize,xmin,xmax);
k660_dms_obs=M;u10ndo=X;sigdo=SD;


figure;semilogy(X,taucxmn,'--o',X,tauibn,'-x',X,taubn,'-d');xlabel('u_{10n} (m/s)');ylabel('<tau>');ylim([0 2.2]);legend('<cov>','ID','C36');grid;title('HIWINGS');ylim([1e-3 3])

figure;
subplot(2,1,1);plot(u10ncu,k660_co2_um,'--o',u10ncw,k660_co2_wm,'-x',u10ndu,k660_co2_obs,'-d');xlabel('u_{10n} (m/s)');ylabel('<k_{660}> (cm/hr)');legend('co2 u','co2 w','co2 data');grid;title('HIWINGS');
subplot(2,1,2);plot(u10ncu,k660_dms_um,'--o',u10ncw,k660_dms_wm,'-x',u10ndu,k660_dms_obs,'-d');xlabel('u_{10n} (m/s)');ylabel('<k_{660}> (cm/hr)');legend('dms u','dms w','dms data');grid


