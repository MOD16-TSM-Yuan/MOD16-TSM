%==============================%
%    Revised Rs from the MOD16 algorithm   %
%==============================%

% Code: Ling Yuan,  Institute of Tibetan Plateau Research, Chinese Academy of Sciences, Beijing, China.
% Questions to: yuanling@itpcas.ac.cn

% References:
% Cleugh et al., 2007. Regional evaporation estimates from flux tower and MODIS satellite data
% Mu et al., 2007. Development of a global evapotranspiration algorithm based on MODIS and global meteorology data
% Mu et al., 2011. Improvements to a MODIS global terrestrial evapotranspiration algorithm
%--------------------------------------------------------------------------

Data=xlsread('E:\Pre_Data\Ali\Test\****.xlsx');
[m,n]=size(Data);

beta   = 2;          % (hPa)
T0     = 273.15;  % zero Kelvin [C]
P0     = 1013.25;  % Standard pressure (hPa)
rd      = 287.0;      % Gas Constant for Dry air, from table 2.1 P25 of Brutsaert 2005 (J kg-1 K-1)
k       = 0.41;       % von Karman constant
g        = 9.8;       % Gravity accelaration (kg s-2)

zh    = Data(:,4);        % reference level of air temperature (m)
zm    = zh;                 
h     = Data(:,5);         % canopy height, (m)
u     = Data(:,6);         % Wind speed, (m/s)
Ta    = Data(:,7);        % air temperature, (C)
Rh    = Data(:,8);        % relative humidity, (%)
Pa    = Data(:,9);        % air pressure, (hPa)
NDVI  = Data(:,10);  % the normalized difference vegetation index
Rn    = Data(:,11);     % net radiation, (W m-2)
G     = Data(:,12);      % soil heat flux, (W m-2)
uf    = Data(:,13);       % frictional velocity
L     = Data(:,14);       % the length of Monin-Obukhov (M-O length)
LE    = Data(:,15);     % LE: latent heat flux, (W m-2)

Cp         = 1013;            % Specific heat (J kg-1 C-1)
rddcp    = rd/Cp;
eps        = 0.622;            % e (unitless) is the ratio of molecular weight of water to dry air (equal to 0.622)
sigma    = 5.6703e-8;   % Stefen-Boltzman's constant, (W m-2 K-4)
lambda  = (2501-2.361*Ta)*1000;                       % Latent heat of water vaporization, (J kg-1)
es          = 6.108*exp((17.27*Ta)./(237.3+Ta));       % Saturation vapour pressure at Ta, (hPa)      
ea          = es.*Rh/100;                                             % Actual vapour pressure (hPa)
VPD      = es-ea;                                                    % VPD: vapor pressure deficit, (hPa)
delta      = (4098 .* es) ./ ((Ta + 237.3).^2);            % Slope of saturation vapour pressure curve at Ta (hPa/degC)
gamma   = (Cp .* Pa) ./ (eps * lambda);              % Psychrometric constant (hPa/degC)
rcorr      = 1.0./((P0./Pa).*(((Ta+273.15)/293.15).^1.75));
rho_a     = 0.3846*Pa./(Ta+273.15);                                          % Density of air (kg m-3)

LAI       = single(real((sqrt(NDVI.*(1+NDVI)./(1-NDVI)))));   % leaf area index
fc           = ((NDVI-0.005)./(0.95-0.005)).^2;                             % vegetation cover fraction                        
nu          = 1.327e-5*(P0./Pa).*(((Ta+T0)/ T0).^1.754);            % kinematic viscousity 

fwet      = (Rh/100).^4;   % fwet: wet surface fraction
fwet(Rh < 70) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% z0h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0m     = 0.138*h;          % aerodynamic roughness length
Tf      = Ta.*(uf.^2)./(k*g*L);    % frictional tempertature
z0h     = (70*nu./uf).*exp(-7.2*((uf.^0.5).*(abs(Tf).^0.25)));  % thermal roughness length
z0h     = min(zh/10,max(z0h,1.0E-10));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m
    if isnan(L(i,1))==0
       if L(i,1)>0
          x(i,1)=0;
          x0(i,1)=0;
          y(i,1)=0;
          y0(i,1)=0;
          Fm(i,1)=-5.3*(zm(i,1)-z0m(i,1))/L(i,1);
          Fh(i,1)=-8.0*(zh(i,1)-z0h(i,1))/L(i,1);
       end
       if L(i,1)==0
          x(i,1)=0;
          x0(i,1)=0;
          y(i,1)=0;
          y0(i,1)=0;
          Fm(i,1)=0;
          Fh(i,1)=0;
       end
       if L(i,1)<0
          x(i,1)=(1-19*zm(i,1)/L(i,1))^0.25;
          x0(i,1)=(1-19*z0m(i,1)/L(i,1))^0.25;
          y(i,1)=(1-11.6*zh(i,1)/L(i,1))^0.5;
          y0(i,1)=(1-11.6*z0h(i,1)/L(i,1))^0.5;
          Fm(i,1)=2*log((1+x(i,1))/(1+x0(i,1)))+log((1+x(i,1)*x(i,1))/(1+x0(i,1)*x0(i,1)))-2*atan(x(i,1))+2*atan(x0(i,1));
          Fh(i,1)=2*log((1+y(i,1))/(1+y0(i,1)));
       end
    else
       x(i,1)=NaN;
       x0(i,1)=NaN;
       y(i,1)=NaN;
       y0(i,1)=NaN;
       Fm(i,1)=NaN;
       Fh(i,1)=NaN;
    end
end
i=0;
Ra  = log(zh./z0h-Fh).*log(zm./z0m-Fm)./(u*k*k);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%
% Radiation partition
Ac     = fc .* Rn;               % the part of Net Radiation allocated to the canopy, (W m-2)
As     = (1 - fc) .* Rn - G; % the part of Net Radiation allocated to soil surface (W m-2)

%%  Rs  %
A   = (delta.*As+Cp*rho_a.*VPD./Ra).*(Rh/100).^(VPD ./ beta);
Rs  = (A-delta).*Ra./gamma-Ra; 



