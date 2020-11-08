% Data=xlsread('D:\Test\Pre_Data\US-DK1\USDK1.xlsx');
Data=xlsread('D:\Test\MODIS\US-Wkg\USWkg.xlsx');
% Data=xlsread('D:\Test\MODIS\US-Fwf\USFwf.xlsx');
% Data=xlsread('D:\Test\MODIS\BJ\BJ.xlsx');
% Data=xlsread('D:\Test\MODIS\Duolun\Duolun.xlsx');
[m,n]=size(Data);

beta=2;
T0=273.15;
P0=1013.25; 
rd=287.0;
h=0.03;
k=0.41;
g=9.8;
ksac=0.622;


gcu=0.00001;
glsh=0.02;
gl_e_wv=0.02;
rblmax=95;
rblmin=60;
Rv=461.5;
VPD_close=36;
VPD_open=6.5;
T_open=12.02;
T_close=-8;

zh=Data(:,4);
zm=zh;
d0=Data(:,5);
U=Data(:,6);
Ta=Data(:,7);
RH=Data(:,8);
P=Data(:,9);
NDVI=Data(:,10);
Rn=Data(:,11);
G0=Data(:,12);
uf=Data(:,13);
H=Data(:,14);
LE=Data(:,15);
SWC=Data(:,16);
FPAR=Data(:,17);
% LAI=Data(:,18);
% VPD=Data(:,19);

Cp=4.2*242; %J*kg-1*K
lamada=(2501-2.361*Ta)*1000; %J*kg-1
rho_a=0.3846*P./(Ta+273.15);
gama=Cp*P./(0.622*lamada);
Es=6.108*exp((17.27*Ta)./(237.3+Ta));%hpa
Ea=Es.*RH/100;%hpa
VPD=Es-Ea;
deta=4098*Es./((Ta+237.3).^2);%hpa/0C
LAI=single(real((sqrt(NDVI.*(1+NDVI)./(1-NDVI))))); 
FPAR=((NDVI-0.05)./(0.95-0.05)).^2;

A=Rn-G0;
Ac=FPAR.*(Rn-G0);
As=(1-FPAR).*(Rn-G0);

for i=1:m
  if isnan(RH(i,1))==0  
    if (RH(i,1)<70)
       Fw(i,1)=0;
    end
    if RH(i,1)>=70
       Fw(i,1)=((RH(i,1)/100).^4);
    end
  end
  if isnan(RH(i,1))==1
     Fw(i,1)=NaN;
  end
end
i=0;
Fs=(RH/100).^(VPD/beta);  

%%%%%%%%%%%%%%%%%%%%%%improved model%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tf=-H./(Cp*rho_a.*uf);
z0m(NDVI>=0.2)=h/8;
z0m(NDVI<0.2)=0.0005;
z0m=z0m';
L=(Ta+T0).*(uf.^2)./(k*g*Tf);
nu=1.327e-5*(P0./P).*(((Ta+T0)/ T0).^1.754); 
z0h=(70*nu./uf).*exp(-7.2*((uf.^0.5).*(abs(Tf).^0.25)));
rcorr=1.0./((P0./P).*(((Ta+273.15)/293.15).^1.75));

for i=1:m
  if isnan(L(i,1))==0
     if L(i,1)>0
        Fm(i,1)=-5.3*(zm(i,1)-z0m(i,1))/L(i,1);
        Fh(i,1)=-8.0*(zh(i,1)-z0h(i,1))/L(i,1);
        x(i,1)=0;
        x0(i,1)=0;
        y(i,1)=0;
        y0(i,1)=0;
      end
      if L(i,1)==0
         Fm(i,1)=0;
         Fh(i,1)=0;
        x(i,1)=0;
        x0(i,1)=0;
        y(i,1)=0;
        y0(i,1)=0;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SWC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% w=0;
% t=0;
% for i=1:48:m
%   w=w+1;
%   SM=swc(48*(w-1)+1:48*w);
%   Sm=min(SM);
%   for t=(48*(w-1)+1):48*w
%      SWC(t,1)=Sm;
%   end
% end
% w=0;
% t=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ra-Rs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ra=log(zh./z0h-Fh).*log(zm./z0m-Fm)./(U*k*k);
Rss=exp(7.28-4.15*(SWC/40.46));%QOMO-USWkg
% Rss=exp(7.87-4.15*(SWC/67));%Arou-USFwf
% Rss=exp(7.99-3.4*(SWC/65));%Maqu-USDK1
% Rss=exp(7.21-16.12*(SWC/44.8));%NAMOR-BJ
% Rss=exp(8.3-10.4*(SWC/39.22));%NASDE-Duolun

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ta-VPD-limited%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:m
     if isnan(Ta(i,1))==0
       if Ta(i,1)>=T_open
           m_Ta(i,1)=1.0;
       end
       if (Ta(i,1)<T_open)&&(Ta(i,1)>T_close)
           m_Ta(i,1)=(Ta(i,1)-T_close)/(T_open-T_close);
       end
       if Ta(i,1)<=T_close
           m_Ta(i,1)=0.1;
       end
     else
       m_Ta(i,1)=NaN;         
     end    
end

for i=1:m
  if isnan(VPD(i,1))==0
    if VPD(i,1)<=VPD_open
       m_VPD(i,1)=1.0;
       Rtotc(i,1)=rblmax;
    end
    if (VPD(i,1)<VPD_close)&&(VPD(i,1)>VPD_open)
       m_VPD(i,1)=(VPD_close-VPD(i,1))/(VPD_close-VPD_open); 
       Rtotc(i,1)=rblmax-((rblmax-rblmin)*(VPD_close-VPD(i,1)))/(VPD_close-VPD_open);
    end
    if (VPD(i,1)>=VPD_close)
       m_VPD(i,1)=0.1;
       Rtotc(i,1)=rblmin;
    end
  else
    m_VPD(i,1)=NaN;
    Rtotc(i,1)=NaN;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ta-VPD-limited%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Yuan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LEs_pot=(deta.*As+Cp*rho_a.*VPD.*(1-FPAR)./Ra).*Fs.*(1-Fw)./(deta+gama.*(1+Rss./Ra));
LEs_wet=(deta.*As+Cp*rho_a.*VPD.*(1-FPAR)./Ra).*Fw./(deta+gama.*(1+Rss./Ra));
% LEs_pot=(deta.*As+Cp*rho_a.*VPD./Ra).*Fs.*(1-Fw)./(deta+gama.*(1+Rss./Ra));
% LEs_wet=(deta.*As+Cp*rho_a.*VPD./Ra).*Fw./(deta+gama.*(1+Rss./Ra));
Gs=0.0038*m_Ta.*m_VPD.*LAI;
Rcs=1./Gs;
LEc=((deta.*Ac+rho_a.*Cp.*VPD./Ra).*FPAR.*(1-Fw))./(deta+gama.*(1+Rcs./Ra));
LEc_wet=(deta.*Ac+Cp*rho_a.*VPD./Ra).*FPAR.*Fw./(deta+gama.*(1+Rcs./Ra));
% LEc=((deta.*Ac+rho_a.*Cp.*VPD./Ra).*(1-Fw))./(deta+gama.*(1+Rcs./Ra));
% LEc_wet=(deta.*Ac+Cp*rho_a.*VPD./Ra).*Fw./(deta+gama.*(1+Rcs./Ra));
LEc(NDVI<0.2)=0;
LEm=LEs_pot+LEs_wet+LEc+LEc_wet;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Yuan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MOD16 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%Evaporation from wet canopy surface%%%%%%%%%%%%%%%%%
for i=1:m
   if Fw(i,1)~=0
      rvc(i,1)=(rho_a(i,1)*Cp)/(gl_e_wv*LAI(i,1)*Fw(i,1));
      rhc(i,1)=1/(glsh*LAI(i,1)*Fw(i,1));
      rrc(i,1)=rho_a(i,1)*Cp/(4*(5.67*10^(-8))*(Ta(i,1)+273.15))^3;
      rhrc(i,1)=(rhc(i,1)*rrc(i,1))/(rhc(i,1)+rrc(i,1));
      LEc_wet_m(i,1)=((deta(i,1)*Ac(i,1)+rho_a(i,1)*Cp*VPD(i,1)/rhrc(i,1))*FPAR(i,1)*Fw(i,1))/(deta(i,1)+gama(i,1)*(1+rvc(i,1)/rhrc(i,1)));
   else
      rvc(i,1)=0;
      rhc(i,1)=0;
      rrc(i,1)=0;
      rhrc(i,1)=0;
      LEc_wet_m(i,1)=0;
   end
end
%%%%%%%%%%%%%%%%Plant transpiration%%%%%%%%%%%%%%%%
%Aerodynamic resistance%
rhc=1.0./glsh;
rrc=rho_a*Cp./(4.0*5.67*10^(-8)*(Ta+273.15).^3);
Rac=(rhc.*rrc./(rhc+rrc));

%Surface conductance to transpiration%
Gs1=0.0065*m_Ta.*m_VPD.*LAI.*rcorr;
Gs1(Rn<0)=0;
Gcu=gcu*rcorr;
Gs2=glsh;
for i=1:m
   if  (LAI(i,1)>0)&&(1-Fw(i,1)>0)
      Cc(i,1)=((Gs2*(Gs(i,1)+Gcu(i,1)))/(Gs2+Gs(i,1)+Gcu(i,1)))*LAI(i,1)*(1-Fw(i,1));
   else
      Cc(i,1)=NaN;
   end
end
Rsc=1./Cc;
LEc_m=((deta.*Ac+rho_a.*Cp.*VPD./Rac).*FPAR.*(1-Fw))./(deta+gama.*(1+Rsc./Rac));
LEc_m(NDVI<0.2)=0;

%%%%%%%%%%%%%%%%Evaporation from soil surface%%%%%%%%%%%%%%%%
Rtot= Rtotc.*rcorr;
rrs=rrc;
rhs=Rtot;
Ras=rrs.*rhs./(rrs+rhs);

LEs_wet_m=(deta.*As+Cp*rho_a.*VPD.*(1-FPAR)./Ras).*Fw./(deta+gama.*Rtot./Ras);
LEs_pot_m=(deta.*As+Cp*rho_a.*VPD.*(1-FPAR)./Ras).*Fs.*(1-Fw)./(deta+gama.*(Rtot./Ras));

LE_mu=LEs_pot_m+LEs_wet_m+LEc_wet_m+LEc_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MOD16 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MOD16-Chang%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zom=0.123*h;
zoh=0.1*zom;
Ra_Thom=log((zm-d0)./zom).*log((zh-d0)./zoh)./(k*k*U);
LEc_Chang=((deta.*Ac+rho_a.*Cp.*VPD./Rac).*FPAR.*(1-Fw))./(deta+gama.*(1+Rsc./Ra_Thom));
LEs_wet_Chang=1.26*Fw.*deta.*As./(deta+gama);
LEs_pot_Chang=1.26*Fs.*(1-Fw).*deta.*As./(deta+gama);
LEc_wet_Chang=(deta.*Ac+Cp*rho_a.*VPD./Ra).*FPAR.*Fw./(deta+gama.*(1+Rsc./Ra_Thom));

LE_Chang=LEc_Chang+LEc_wet_Chang+LEs_wet_Chang+LEs_pot_Chang;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MOD16-Chang%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%P-T%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ft=exp(-((Ta-25)/25).^2);
LEV=1.26*(1-Fw).*FPAR.*Ft.*deta.*Ac./(deta+gama);
LEVw=1.26*Fw.*deta.*Ac./(deta+gama);
LES=1.26*(1-Fw).*Fs.*deta.*As./(deta+gama);
LESw=1.26*Fw.*deta.*As./(deta+gama);
LE_PT=LEV+LEVw+LES+LESw;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%P-T%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Daily%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum_LEm=0;
sum_LE=0;
sum_LE_mu=0;
sum_LE_Chang=0;
sum_LE_PT=0;
k1=0;
k2=0;
k3=0;
k4=0;
k5=0;
t=0;

for i=1:48:m
   t=t+1;
   for j=48*(t-1)+1:48*t
       if isnan(LE(j,1))==0
          sum_LE=sum_LE+LE(j,1);
          k1=k1+1;
       end
       if isnan(LEm(j,1))==0
          sum_LEm=sum_LEm+LEm(j,1);
          k2=k2+1;
       end
       if isnan(LE_mu(j,1))==0
          sum_LE_mu=sum_LE_mu+LE_mu(j,1);
          k3=k3+1;
       end
       if isnan(LE_Chang(j,1))==0
          sum_LE_Chang=sum_LE_Chang+LE_Chang(j,1);
          k4=k4+1;
       end
       if isnan(LE_PT(j,1))==0
          sum_LE_PT=sum_LE_PT+LE_PT(j,1);
          k5=k5+1;
       end
   end
    Daily(t,1)=Data(i,1);
    Daily(t,2)=Data(i,2);
    Daily(t,3)=Data(i,3);
    Daily(t,4)=sum_LE/k1;
    Daily(t,5)=k1;
    Daily(t,6)=sum_LEm/k2;
    Daily(t,7)=k2;
    Daily(t,8)=sum_LE_mu/k3;
    Daily(t,9)=k3;
    Daily(t,10)=sum_LE_Chang/k4;
    Daily(t,11)=k4;
    Daily(t,12)=sum_LE_PT/k5;
    Daily(t,13)=k5;
    sum_LEm=0;
    sum_LE=0;
    sum_LE_mu=0;
    sum_LE_Chang=0;
    sum_LE_PT=0;
    k1=0;
    k2=0;
    k3=0;
    k4=0;
    k5=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Daily%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Average-Month%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=[31,28,31,30,31,30,31,31,30,31,30,31];
sum=0;
k0=0;
sum_m=0;
kk=0;
sum_mu=0;
km=0;
sum_Chang=0;
kc=0;
sum_PT=0;
kp=0;
day=0;
for f=1:12
 for i=1:48
  for j=1:Y(f)
       if isnan(Data(48*(day+j-1)+i,15))==0
          sum=sum+Data(48*(day+j-1)+i,15); 
          k0=k0+1;
       end
       if isnan(LEm(48*(day+j-1)+i,1))==0
          sum_m=sum_m+LEm(48*(day+j-1)+i,1); 
          kk=kk+1;
       end
       if isnan(LE_mu(48*(day+j-1)+i,1))==0
          sum_mu=sum_mu+LE_mu(48*(day+j-1)+i,1); 
          km=km+1;
       end
       if isnan(LE_Chang(48*(day+j-1)+i,1))==0
          sum_Chang=sum_Chang+LE_Chang(48*(day+j-1)+i,1); 
          kc=kc+1;
       end
       if isnan(LE_PT(48*(day+j-1)+i,1))==0
          sum_PT=sum_PT+LE_PT(48*(day+j-1)+i,1); 
          kp=kp+1;
       end
  end
  Ay(48*(f-1)+i,1)=sum/k0;
  Ay(48*(f-1)+i,2)=k0;
  Ay(48*(f-1)+i,3)=sum_m/kk;
  Ay(48*(f-1)+i,4)=kk;
  Ay(48*(f-1)+i,5)=sum_mu/km;
  Ay(48*(f-1)+i,6)=km;
  Ay(48*(f-1)+i,7)=sum_Chang/kc;
  Ay(48*(f-1)+i,8)=kc;
  Ay(48*(f-1)+i,9)=sum_PT/kp;
  Ay(48*(f-1)+i,10)=kp;
  sum=0;
  k0=0;
  sum_m=0;
  kk=0;
  sum_mu=0;
  km=0;
  sum_Chang=0;
  kc=0;
  sum_PT=0;
  kp=0;
 end
 day=day+Y(f);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Average-Month%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Tongji-Half-hour%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1=0;
z2=0;
z3=0;
z4=0;
SMB1=0;
SMB2=0;
SMB3=0;
SMB4=0;
for i=1:m
    if (isnan(LE(i,1))==0)&&(isnan(LEm(i,1))==0)
        z1=z1+1;
        Y_LE1(z1,1)=Data(i,1);
        Y_LE1(z1,2)=Data(i,2);
        Y_LE1(z1,3)=Data(i,3);
        Y_LE1(z1,4)=Data(i,19);
        Y_LE1(z1,5)=Data(i,20);
        Y_LE1(z1,6)=LE(i,1);
        Y_LEm(z1,1)=LEm(i,1);
    end
    if (isnan(LE(i,1))==0)&&(isnan(LE_mu(i,1))==0)
        z2=z2+1;
        Y_LE2(z2,1)=LE(i,1);
        Y_LE_mu(z2,1)=LE_mu(i,1);
    end
    if (isnan(LE(i,1))==0)&&(isnan(LE_Chang(i,1))==0)
        z3=z3+1;
        Y_LE3(z3,1)=LE(i,1);
        Y_LE_Chang(z3,1)=LE_Chang(i,1);
    end
    if (isnan(LE(i,1))==0)&&(isnan(LE_PT(i,1))==0)
        z4=z4+1;
        Y_LE4(z4,1)=LE(i,1);
        Y_LE_PT(z4,1)=LE_PT(i,1);
    end
end
%%%%%%%R2%%%%%%%%%%%%
C1=corrcoef(Y_LE1(:,6),Y_LEm);
C2=corrcoef(Y_LE2,Y_LE_mu);
C3=corrcoef(Y_LE3,Y_LE_Chang);
C4=corrcoef(Y_LE4,Y_LE_PT);
%%%%%%%xielv%%%%%%%%%%%%
p1=polyfit(Y_LE1(:,6),Y_LEm,1);
p2=polyfit(Y_LE2,Y_LE_mu,1);
p3=polyfit(Y_LE3,Y_LE_Chang,1);
p4=polyfit(Y_LE4,Y_LE_PT,1);
for jj=1:z1
    SMB1=SMB1+Y_LE1(jj,6)-Y_LEm(jj,1);
    MB1(jj,1)=Y_LE1(jj,6)-Y_LEm(jj,1);
end
for jj=1:z2
    SMB2=SMB2+Y_LE2(jj,1)-Y_LE_mu(jj,1);
    MB2(jj,1)=Y_LE2(jj,1)-Y_LE_mu(jj,1);
end
for jj=1:z3
    SMB3=SMB3+Y_LE3(jj,1)-Y_LE_Chang(jj,1);
    MB3(jj,1)=Y_LE3(jj,1)-Y_LE_Chang(jj,1);
end
for jj=1:z4
     SMB4=SMB4+Y_LE4(jj,1)-Y_LE_PT(jj,1);
     MB4(jj,1)=Y_LE4(jj,1)-Y_LE_PT(jj,1);
end
T(1,1)=C1(1,2);
T(2,1)=C2(1,2);
T(3,1)=C3(1,2);
T(4,1)=C4(1,2);
T(1,2)=p1(1);
T(2,2)=p2(1);
T(3,2)=p3(1);
T(4,2)=p4(1);
T(1,3)=p1(2);
T(2,3)=p2(2);
T(3,3)=p3(2);
T(4,3)=p4(2);
T(1,4)=sqrt(mean((Y_LE1(:,6)-Y_LEm).^2));
T(2,4)=sqrt(mean((Y_LE2-Y_LE_mu).^2));
T(3,4)=sqrt(mean((Y_LE3-Y_LE_Chang).^2));
T(4,4)=sqrt(mean((Y_LE4-Y_LE_PT).^2));
T(1,5)=SMB1/z1;
T(2,5)=SMB2/z2;
T(3,5)=SMB3/z3;
T(4,5)=SMB4/z4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Tongji-Day%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y1=0;
y2=0;
y3=0;
MD1=0;
MD2=0;
MD3=0;
for i=1:365
    if (isnan(Daily(i,4))==0)&&(isnan(Daily(i,6))==0)&&Daily(i,5)>40&&Daily(i,7)>40
        y1=y1+1;
        D_LE1(y1,1)=Daily(i,1);
        D_LE1(y1,2)=Daily(i,2);
        D_LE1(y1,3)=Daily(i,3);
        D_LE1(y1,4)=Daily(i,4);
        D_LEm(y1,1)=Daily(i,6);
    end
    if (isnan(Daily(i,4))==0)&&(isnan(Daily(i,8))==0)
        y2=y2+1;
        D_LE2(y2,1)=Daily(i,4);
        D_LE_mu(y2,1)=Daily(i,8);
    end
    if (isnan(Daily(i,4))==0)&&(isnan(Daily(i,8))==0)
        y3=y3+1;
        D_LE3(y3,1)=Daily(i,4);
        D_LE_PT(y3,1)=Daily(i,12);
    end
end
%%%%%%%R2%%%%%%%%%%%%
CD1=corrcoef(D_LE1(:,4),D_LEm);
CD2=corrcoef(D_LE2,D_LE_mu);
CD3=corrcoef(D_LE3,D_LE_PT);
%%%%%%%xielv%%%%%%%%%%%%
pd1=polyfit(D_LE1(:,4),D_LEm,1);
pd2=polyfit(D_LE2,D_LE_mu,1);
pd3=polyfit(D_LE3,D_LE_PT,1);

G1=mean(D_LE1(:,4));
G2=mean(D_LE2);
G3=mean(D_LE3);
AMm=mean(D_LEm);
AM_mu=mean(D_LE_mu);
AM_PT=mean(D_LE_PT);
fz1=0;
fz2=0;
fz3=0;
fm1=0;
fm2=0;
fm3=0;

for jj=1:y1
    MD1=MD1+D_LE1(jj,4)-D_LEm(jj,1);
    fz1=fz1+(D_LEm(jj,1)-D_LE1(jj,4))^2;
    fm1=fm1+(abs(D_LEm(jj,1)-G1)+abs(D_LE1(jj,4)-G1))^2;
end
for jj=1:y2
    MD2=MD2+D_LE2(jj,1)-D_LE_mu(jj,1);
    fz2=fz2+(D_LE_mu(jj,1)-D_LE2(jj,1))^2;
    fm2=fm2+(abs(D_LE_mu(jj,1)-G2)+abs(D_LE2(jj,1)-G2))^2;
end
for jj=1:y3
    MD3=MD3+D_LE3(jj,1)-D_LE_PT(jj,1);
    fz3=fz3+(D_LE_PT(jj,1)-D_LE3(jj,1))^2;
    fm3=fm3+(abs(D_LE_PT(jj,1)-G3)+abs(D_LE3(jj,1)-G3))^2;
end
TD(1,1)=CD1(1,2);
TD(2,1)=CD2(1,2);
TD(3,1)=CD3(1,2);
TD(1,2)=pd1(1);
TD(2,2)=pd2(1);
TD(3,2)=pd3(1);
TD(1,3)=pd1(2);
TD(2,3)=pd2(2);
TD(3,3)=pd3(2);
TD(1,4)=sqrt(mean((D_LE1(:,4)-D_LEm).^2));
TD(2,4)=sqrt(mean((D_LE2-D_LE_mu).^2));
TD(3,4)=sqrt(mean((D_LE3-D_LE_PT).^2));
TD(1,5)=MD1/y1;
TD(2,5)=MD2/y2;
TD(3,5)=MD3/y3;
TD(1,6)=1-fz1/fm1;
TD(2,6)=1-fz2/fm2;
TD(3,6)=1-fz3/fm3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Tongji-Month%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=0;
x2=0;
x3=0;
Mx1=0;
Mx2=0;
Mx3=0;
for i=1:365
    if (isnan(Ay(i,1))==0)&&(isnan(Ay(i,3))==0)&&Ay(i,2)>20&&Ay(i,4)>20
        x1=x1+1;
        M_LE1(x1,1)=Ay(i,1);
        M_LE_m(x1,1)=Ay(i,3);
    end
    if (isnan(Ay(i,1))==0)&&(isnan(Ay(i,5))==0)
        x2=x2+1;
        M_LE2(x2,1)=Ay(i,1);
        M_LE_mu(x2,1)=Ay(i,5);
    end
    if (isnan(Ay(i,1))==0)&&(isnan(Ay(i,9))==0)
        x3=x3+1;
        M_LE3(x3,1)=Ay(i,1);
        M_LE_PT(x3,1)=Ay(i,9);
    end
end
%%%%%%%R2%%%%%%%%%%%%
CM1=corrcoef(M_LE1,M_LE_m);
CM2=corrcoef(M_LE2,M_LE_mu);
CM3=corrcoef(M_LE3,M_LE_PT);
%%%%%%%xielv%%%%%%%%%%%%
pm1=polyfit(M_LE1,M_LE_m,1);
pm2=polyfit(M_LE2,M_LE_mu,1);
pm3=polyfit(M_LE3,M_LE_PT,1);

GM1=mean(M_LE1);
GM2=mean(M_LE2);
GM3=mean(M_LE3);
Mm=mean(M_LE_m);
M_mu=mean(M_LE_mu);
M_PT=mean(M_LE_PT);
Mz1=0;
Mz2=0;
Mz3=0;
Mm1=0;
Mm2=0;
Mm3=0;


for jj=1:x1
    Mx1=Mx1+M_LE1(jj,1)-M_LE_m(jj,1);
    Mz1=Mz1+(M_LE_m(jj,1)-M_LE1(jj,1))^2;
    Mm1=Mm1+(abs(M_LE_m(jj,1)-GM1)+abs(M_LE1(jj,1)-GM1))^2;
end
for jj=1:x2
    Mx2=Mx2+M_LE2(jj,1)-M_LE_mu(jj,1);
    Mz2=Mz2+(M_LE_mu(jj,1)-M_LE2(jj,1))^2;
    Mm2=Mm2+(abs(M_LE_mu(jj,1)-GM2)+abs(M_LE2(jj,1)-GM2))^2;
end
for jj=1:x2
    Mx3=Mx3+M_LE3(jj,1)-M_LE_PT(jj,1);
    Mz3=Mz3+(M_LE_PT(jj,1)-M_LE3(jj,1))^2;
    Mm3=Mm3+(abs(M_LE_PT(jj,1)-GM3)+abs(M_LE3(jj,1)-GM3))^2;
end
Tm(1,1)=CM1(1,2);
Tm(2,1)=CM2(1,2);
Tm(3,1)=C3(1,2);
Tm(1,2)=pm1(1);
Tm(2,2)=pm2(1);
Tm(3,2)=pm3(1);
Tm(1,3)=pm1(2);
Tm(2,3)=pm2(2);
Tm(3,3)=pm3(2);
Tm(1,4)=sqrt(mean((M_LE1-M_LE_m).^2));
Tm(2,4)=sqrt(mean((M_LE2-M_LE_mu).^2));
Tm(3,4)=sqrt(mean((M_LE3-M_LE_PT).^2));
Tm(1,5)=Mx1/x1;
Tm(2,5)=Mx2/x2;
Tm(3,5)=Mx3/x3;
Tm(1,6)=1-Mz1/Mm1;
Tm(2,6)=1-Mz2/Mm2;
Tm(3,6)=1-Mz3/Mm3;


