% clear all
% close all
cd c:\data\cwf\dropbox\matlabstf\cwf\whoi_bulk;
% Load in the GASEX 10 minute means

load AveMet10new;%load in data set from Gasex2008

% Load in Wades delta CO2 time series
% 
% delta CO2 units: mmol/m^3

load DeltaCO2        % yday10 dCO210
%****   deltaC defined as Xw-alpha*Xa   so Flux = k*deltaC    *************

% Run vectorized version of COAREG
ydms=load('SO_GASEX_UH_DMS_Flux_Hourly_Ver2b_nohds2.txt');
%Time_GMT	FluxDuration_sec	DMS_pptv_air	DMS_uM_m3_air	DMS_uM_m3_sw	DMSflux_uM_m2_d	DMSflux_error	kDMS_cm_hr	kDMS_error	U10N_m_s	Ustar_COARE_m_s	RWspd_m_s	Sc_DMS	RWdir_deg	Tair_C	SST_C	Sal_ppth	Lat	Lon	SOG_kts	Gyro_deg	Longwave_W_m2	Shortwave_W_m2	z_L
jd_dms_H=ydms(:,1)-ydms(1,1)+62.5;
zdms=interp1(jd_dms_H,ydms,yday10);
dDMS=zdms(:,5)/1e3;%DMS water concentration in millimole/m^3
u=U1010;t=T1010;rh=RH1010;P=Pair10;ts=SST10;Rs=Solardn10;Rl=IRdn10;lat=Lat10;zi=700;rain=0;Ss=35;cp=NaN*ones(1,length(u));sigH=cp;zrf_u=10;zrf_t=10;zrf_q=10;
zu=10*ones(1,length(U1010))';zt=10*ones(1,length(U1010))';zq=10*ones(1,length(U1010))';  %height of wind speed, air temperature, relative humdity data
[A G] = coareG36vn_co2(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q,dCO210);%)CO2(U1010,zu,T1010,zt,RH1010,zh,Pair10,SST10,Solardn10,IRdn10,Lat10,700,dCO210);
[Ad Gd] = coareG36vn_dms(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q,dDMS);
[Ao Go] = coareG36vn_ozo(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q,dDMS,cp);
[Ahe Ghe] = coareG36vn_hel(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q,0*dDMS);
[Asf Gsf] = coareG36vn_sf6(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q,0*dDMS);

% The COAREG algorithm gives k in m/s
% del_CO2 is given in mmol/m^3
% therefore, the flux below is given in mmol/m^2/s 
% multiplication by 44 grams/mole is not necessary
%
% wCO2 => m/s mmol/m^3 => mmol/m2/s
%
% Note that the fluxes from Gasex08 are in units of:
%
% wCO2 => mg/m2/s
%
% we then divide by the molecular mass to get:
%
% wCO2/molecular mass => mg/m2/s mol/gram => mmol/m2/s

wCO2 = G(:,12);
wDMS=Gd(:,12);
% 

figure;plot(yday10,wDMS*365,'b',yday10,zdms(:,6)*365/1e3,'or')
%axis([65 101 -35 5])
xlabel('Yearday 2008','fontsize',16)
ylabel('F_{DMS} [mmol/m^2/yr]','fontsize',16);
legend('COARE36G DMS','Eddy Covariance');

% figure;plot(G(:,13),Gh(:,13),'.',G(:,13),Gs(:,13),'.',[0 8e-4],[0 8e-4],'--');
% xlabel('kCO2_{660} vector [m/s]','fontsize',16)
% ylabel('k_{660}  [m/s]','fontsize',16);
% legend('COARE30G_{He}','COARE30G_{SF6}');

figure;semilogy(yday10,G(:,5)*1e2*3600,yday10,Gd(:,5)*1e2*3600,'.');
xlabel('Yearday 2008','fontsize',16)
ylabel('k  [cm/h]','fontsize',16);
legend('COARE36G_{CO2}','COARE36G_{DMS}');

figure;semilogy(U1010,G(:,14)*1e2*3600,'.',U1010,Gd(:,14)*1e2*3600,'.',U1010,Ghe(:,13)*1e2*3600,'.',U1010,Gsf(:,13)*1e2*3600,'.');axis([0 25 1 500]);
%figure;plot(U1010,G(:,14)*1e2*3600,'.',U1010,Gd(:,14)*1e2*3600,'.',U1010,Ghe(:,13)*1e2*3600,'.',U1010,Gsf(:,13)*1e2*3600,'.');axis([0 25 1 500]);
ylabel('k_{660}  [cm/h]','fontsize',16);xlabel('U10 (m/s)','fontsize',16);grid
legend('COARE36G_{CO2}','COARE36G_{DMS}','COARE36G_{3He}','COARE36G_{SF6}');
