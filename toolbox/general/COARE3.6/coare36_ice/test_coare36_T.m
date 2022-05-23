%%%%%%%%%  Program to exercise the coare3.6
%allows water salinity as an input
%if sst<freezing point assume surface is ice.
%This version sets sigH as a simple quadratic in windspeed and uses a fixed
%wave age.

%%%%%%%%%%%%   Set bulk parameters  %%%%%%%%%%%
jj=2:20;
n=length(jj);
u=5*ones(1,n);%wind speed           %wind speed 2 to 20 m/s
zu=18*ones(1,n);%ht of wiind speed
t=0*ones(1,n);%air temp (C)
zt=15*ones(1,n);%height of air tem
rh=80*ones(1,n);%relative humidity
zq=15*ones(1,n);%height of RH
P=1010*ones(1,n);%pressure (mb)
ts=-5+.5*jj;%surface temperature, -5 to +5 C
Rs=100*ones(1,n);%downward solar flux
Rl=400*ones(1,n);%downward IR flux
lat=20*ones(1,n);%latitude
zi=600*ones(1,n);%mixed layer height
rain=0*ones(1,n);%rain rate (mm/hr)
Ss=15*ones(1,n);%salinity (PSU)

%%%%%%%%  set height for 2nd means 
zrf_u=5;
zrf_t=5;
zrf_q=5;
%%%%%%%%%%%%%%
%%%%%%%%%%%  specify sig wave height and wave age
sigH=max(.25,.015*u.^2);%significant wave ht
Wa=1.3;%wave age (Cp/U10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp=Wa*u;
hsig=(0.02*(cp./u).^1.1-0.0025).*u.^2;%significant wave ht wave model parameterization
hsig=max(hsig,.25);
%sigH=hsig;
Tf=-0.0575*Ss+1.71052E-3*Ss.^1.5-2.154996E-4*Ss.*Ss;%water freezing point (C)
%%%%%%%%   call bulk algo %%%%%%    if ts<Tf algorithm assumes ice
A=coare36vn_zrf(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,Ss,cp,sigH,zrf_u,zrf_t,zrf_q);
%A=[usrw tauw hsbw hlbw hbbw hsbbw hlwebbw tsrw qsrw zow  zotw zoqw Cdw Chw Cew  Lw zetw dter dqer tkt Urf Trf Qrf RHrf UrfN Rnl Le rhoa UN U10 U10N Cdn_10 Chn_10 Cen_10 RF Qs Evap T10 Q10 RH10 ug Whf Edis];
%    1     2   3     4   5     6     7      8     9   10   11  12   13  14   15  16  17   18   19  20  21  22  23  24   25   26  27  28  29  30  31   32     33   34     35 36  37   38  39  40  41  42   43 
usrw=A(:,1);% u*
taub=A(:,2);%stress
hsb=A(:,3);%sensible heat
hlb=A(:,4);%latent heat
zo=A(:,10);%velocity roughness
cd=A(:,13);%ambient drag coeff at measure ht
ch=A(:,14);%ambient heat xfer coeff at measure ht
dter=A(:,18);%cool skin
RHrf=A(:,24);RH10=A(:,40);%Rh at 2nd height and 10m
u10n=A(:,31);%10m neutral wind
cdn_10=A(:,32);%neutral 10m drag coeff
chn_10=A(:,33);%neutral 10m heat xfer coeff
Whf=A(:,42);%white cap fraction
Edis=A(:,43);%wave breaking energy dissipation

B=coare35vn(u,zu,t,zt,rh,zq,P,ts,Rs,Rl,lat,zi,rain,cp,sigH);


 figure;
 subplot(3,1,1);plot(ts,taub);xlabel('ts (C)');ylabel('Surface stress (N/m^2)');grid;title(['Salinity= ' num2str(Ss(1)) ' U= ' num2str(u(1)) ' Ta= ' num2str(t(1)) ' RH= ' num2str(rh(1)) ' Wave Age= ' num2str(Wa(1))]); 
 subplot(3,1,2);plot(ts,hsb);xlabel('ts (C)');ylabel('Sensible Heat Flux (W/m^2)');grid;
 subplot(3,1,3);plot(ts,hlb);xlabel('ts (C)');ylabel('Latent Heat Flux (W/m^2)');grid;

%  figure;
%  subplot(2,1,1);semilogy(u10n,Whf*100);xlabel('U10n (m/s)');ylabel('Whitecap Fraction (%)');grid;;ylim([.1 10]);
%  title(['Ts= ' num2str(ts(1)) ' Ta= ' num2str(t(1)) ' RH= ' num2str(rh(1)) ' Wave Age= ' num2str(Wa(1))]); 
%  subplot(2,1,2);semilogy(u10n,Edis);xlabel('U10n (m/s)');ylabel('Wave Energy Dissipation W/m^2');grid;;ylim([1e-3 10]);
% 
