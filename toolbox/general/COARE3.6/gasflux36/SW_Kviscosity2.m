
function nu = SW_Kviscosity2(T,S)
% from SeawaterLab v 3.1.2 07/08/2016 Copyright (C) Kishor G Nayar et al., 2016
% T = SST, C
% S = salinity ppth
% output: kinematic viscosity, m^2/s

P0 = SW_Psat2(T,S)/1E6;
P0(T<100) = 0.101325;
mu  = SW_Viscosity2(T,S);
rho = SW_Density2(T,S,P0);
nu  = mu./rho;
end
function mu = SW_Viscosity2(T,S)
% from SeawaterLab v 3.1.2 07/08/2016 Copyright (C) Kishor G Nayar et al., 2016
% T = SST, C
% S = salinity ppth
% output, dynamic viscosity, kg/s/m

S = S/1000;
a = [1.5700386464E-01; 6.4992620050E+01; -9.1296496657E+01; 4.2844324477E-05;...
     1.5409136040E+00; 1.9981117208E-02; -9.5203865864E-05; 7.9739318223E+00;...
    -7.5614568881E-02; 4.7237011074E-04];
mu_w = a(4) + 1./(a(1)*(T+a(2)).^2+a(3));
A  = a(5) + a(6) * T + a(7) * T.^2;
B  = a(8) + a(9) * T + a(10)* T.^2;
mu = mu_w.*(1 + A.*S + B.*S.^2);
end
%------------------------------------------------------------------------------------------------------
function rho = SW_Density2(T,S,P)
% from SeawaterLab v 3.1.2 07/08/2016 Copyright (C) Kishor G Nayar et al., 2016
% T = SST, C
% S = salinity ppth
% P = pressure in MPa (Pa/1e6)
% output: kg/m^3

P0 = SW_Psat2(T,S)/1E6;
P0(T<100) = 0.101325;
s = S/1000;
a = [9.9992293295E+02; 2.0341179217E-02; -6.1624591598E-03; 2.2614664708E-05; -4.6570659168E-08];
b = [8.0200240891E+02; -2.0005183488E+00; 1.6771024982E-02; -3.0600536746E-05; -1.6132224742E-05];
rho_w = a(1) + a(2)*T + a(3)*T.^2 + a(4)*T.^3 + a(5)*T.^4;
D_rho = b(1)*s + b(2)*s.*T + b(3)*s.*T.^2 + b(4)*s.*T.^3 + b(5)*s.^2.*T.^2;
rho_sw_sharq   = rho_w + D_rho;

c = [5.0792E-04; -3.4168E-06; 5.6931E-08; -3.7263E-10; 1.4465E-12; -1.7058E-15;...
    -1.3389E-06; 4.8603E-09; -6.8039E-13];
d=[-1.1077e-06; 5.5584e-09; -4.2539e-11; 8.3702e-09];
F_P = exp( (P-P0).*(c(1) + c(2)*T + c(3)*T.^2 + c(4)*T.^3 + c(5)*T.^4 + c(6)*T.^5 + S.*(d(1) + ...
    d(2)*T + d(3)*T.^2)) + 0.5*(P.^2-P0.^2).*(c(7) + c(8)*T + c(9)*T.^3   + d(4)*S));
rho = rho_sw_sharq.*F_P;
end

function Pv = SW_Psat2(T,S) % alternate seawater saturation vapor pressure function
% from SeawaterLab v 3.1.2 07/08/2016 Copyright (C) Kishor G Nayar et al., 2016
% T = SST, C
% S = salinity ppth
% output, saturation water vapor pressure, Pa

T = T + 273.15;
a = [-5.8002206E+03; 1.3914993E+00; -4.8640239E-02; 4.1764768E-05; -1.4452093E-08; 6.5459673E+00];
Pv_w = exp((a(1)./T) + a(2) + a(3)*T + a(4)*T.^2 + a(5)*T.^3 + a(6)*log(T));
b  = [-4.5818 * 10 ^ -4; -2.0443 * 10 ^ -6];
Pv   = Pv_w.*exp(b(1)*S+b(2)*S.^2);
end