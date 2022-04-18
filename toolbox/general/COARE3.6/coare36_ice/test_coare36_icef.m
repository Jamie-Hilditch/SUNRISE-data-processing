%%%%%%%%%  Program to exercise the coare3.6
%allows water salinity as an input
%if sst<freezing point assume surface is ice.

%%%%%%%%%%%%   Set bulk parameters  %%%%%%%%%%%
jj=2:20;
n=length(jj);
u=5*ones(1,n);%wind speed
zu=18*ones(1,n);%ht of wiind speed
t=-4*ones(1,n);%air temp (C)
zt=15*ones(1,n);%height of air tem
rh=98*ones(1,n);%relative humidity
zq=15*ones(1,n);%height of RH
P=1010*ones(1,n);%pressure (mb)
ts=-7+.5*jj;%-5*ones(1,n);%surface temperature
Rs=100*ones(1,n);%downward solar flux
Rl=300*ones(1,n);%downward IR flux
lat=10*ones(1,n);%latitude
zi=600*ones(1,n);%mixed layer height
rain=0*ones(1,n);%rain rate (mm/hr)
Ss=35*ones(1,n);%salinity (PSU)

%%%%%%%%  set height for 2nd means 
zrf_u=5;
zrf_t=5;
zrf_q=5;
%%%%%%%%%%%%%%
sigH=max(.25,.015*u.^2);%significant wave ht
Wa=0.7;%wave age (Cp/U10)
cp=Wa*u;
hsig=(0.02*(cp./u).^1.1-0.0025).*u.^2;%significant wave ht wave model parameterization
hsig=max(hsig,.25);
%sigH=hsig;
Tf=-0.0575*Ss+1.71052E-3*Ss.^1.5-2.154996E-4*Ss.*Ss;%water freezing point (C)

%%%%%%%%   call bulk algo
%%%%%%    if ts<Tf assumes ice
A=coare36vn_ice_zrf_f(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q);
%A=[usrw tauw hsbw hlbw hbbw hsbbw hlwebbw tsrw qsrw zow  zotw zoqw Cdw Chw Cew  Lw zetw dter dqer tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis usri taui hsbi hlbi hbbi hsbbi zeti Cdn_10i Chn_10i Cen_10i Cdi Chi Cei];
%    1     2   3     4   5     6     7      8     9   10   11  12   13  14   15  16  17   18   19  20  21  22  23  24   25   26  27  28  29  30  31   32     33   34     35 36  37   38  39  40  41  42   43  44   45   46   47   48    49   50      51     52  53
usrw=A(:,1);% u*
tauw=A(:,2);
hsbw=A(:,3);%sensible heat
hlbw=A(:,4);%latent heat
zow=A(:,10);
cdw=A(:,13);
chw=A(:,14);
dter=A(:,18);%cool skin
RHrfw=A(:,24);RH10w=A(:,40);%Rh at 2nd height and 10m
cdn_10w=A(:,32);%neutral 10m drag coeff
chn_10w=A(:,33);%neutral 10m heat xfer coeff
u10nw=A(:,31);%10-m neutral wind speed
Whf=A(:,42);%white cap fraction
Edis=A(:,43);%wave breaking energy dissipation
usri=A(:,44);% u*
taui=A(:,45);
hsbi=A(:,46);%sensible heat
hlbi=A(:,47);%latent heat
zeti=A(:,50);
cdn_10i=A(:,51);%neutral 10m drag coeff
chn_10i=A(:,52);%neutral 10m heat xfer coeff
cdi=A(:,54);%neutral 10m drag coeff
chi=A(:,55);%neutral 10m heat xfer coeff



% figure;plot(u10n,Whf);xlabel('U10n (m/s)');ylabel('Whitecap Fraction');grid;
% figure;loglog(u10n,Edis);xlabel('U10n (m/s)');ylabel('Wave Energy Dissipation W/m^2');grid;
figure;plot(ts,hlbw,ts,hsbw,ts,hlbi,'.',ts,hsbi,'.');xlabel('ts (C)');ylabel('Hl & Hs (W/m^2)');grid;title(['U= ' num2str(u(1)) '  Ta= ' num2str(t(1)) '  RH= ' num2str(rh(1)) ' Salinity= ' num2str(Ss(1))]);
figure;plot(ts,cdn_10w,ts,chn_10w,ts,cdn_10i,'o',ts,chn_10i,'.');xlabel('ts (C)');ylabel('Cdn_{10} and Chn_{10}');grid;ylim([1 3]);

hf=.01;%ice free board
f=0.5*ones(n,1);%ice fraction
ii=find(ts>Tf);f(ii)=0;
Dmin=8;Dmax=300;ce=0.3;beta=1;s=05;
Ast=1./(1-(Dmin./Dmax).^(1./beta));
Di=Dmin*(Ast./(Ast-f)).^beta;
Dw=Di.*(1-sqrt(f))./sqrt(f);
Sc=(1-exp(-s*Dw./hf));
zowf=2e-4;%water roughness between ice
zoif=5e-4;
Cdnfloe=f.*hf./Di.*Sc.^2*ce/2.*(log(hf./zowf)).^2./(log(zu'./zowf)).^2;
taufloe=Cdnfloe.*u'.^2./(1-psiu_26(zeti)./log(zu'./zoif)).^2;
%taufloe=taufloe';

figure;
subplot(3,1,1);plot(ts,hsbw.*(1-f)+f.*hsbi,'-o',ts,hsbw,ts,hsbi);xlabel('ts (C)');ylabel('Hs (W/m^2)');grid;title(['U= ' num2str(u(1)) '  Ta= ' num2str(t(1)) '  RH= ' num2str(rh(1)) ' Salinity= ' num2str(Ss(1)) '  f_{ice}= ' num2str(f(1)) ' Freeboard = ' num2str(hf)]);legend('hsb mosaic','hsbw','hsbi');
subplot(3,1,2);plot(ts,hlbw.*(1-f)+f.*hlbi,'-o',ts,hlbw,ts,hlbi);xlabel('ts (C)');ylabel('Hl (W/m^2)');grid;legend('hlb mosaic','hlbw','hlbi');
subplot(3,1,3);plot(ts,tauw.*(1-f)+f.*taui+taufloe,'-o',ts,tauw,ts,taui);xlabel('ts (C)');ylabel('Tau (N/m^2)');grid;legend('tau Combined','taubw','taubi');

