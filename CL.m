% Data=xlsread('E:\Pre_Data\Maqu\Test\Maqu.xlsx');
Data=xlsread('F:\Test\Pre_Data\Maqu\Test\Data.xlsx');
[m,n]=size(Data);

beta=50;
T0=273.15;
P0=1013.25; 
h=0.03;
k=0.41;
g=9.8;
ksac=0.622;
VPD_close=36;
VPD_open=6.5;
T_open=12.02;
T_close=-8;

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
uf=Data(:,10);
LE=Data(:,11);
L=Data(:,12);
Rn=Data(:,13);
G0=Data(:,14);
NDVI=Data(:,17);
FPAR=Data(:,19);
swc=Data(:,15);
LAI=Data(:,18);
% Gs=Data(:,16);
% FPAR(Rn<0)=0;

Cp=4.2*242; %J*kg-1*K
lamada=(2501-2.361*Ta)*1000; %J*kg-1
rho_a=0.3846*P./(Ta+273.15);
gama=Cp*P./(0.622*lamada);
Es=6.108*exp((17.27*Ta)./(237.3+Ta));%hpa
Ea=Es.*RH/100;%hpa
VPD=Es-Ea;
deta=4098*Es./((Ta+237.3).^2);%hpa/0C
nu=1.327e-5*(P0./P).*(((Ta+T0)/ T0).^1.754); 
% LAI=single(real((sqrt(NDVI.*(1+NDVI)./(1-NDVI))))); 
% FPAR=((NDVI-0.05)./(0.95-0.05)).^2;

A=Rn-G0;
As=(1-FPAR).*(Rn-G0);
Ac=FPAR.*(Rn-G0);


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

z0m(NDVI>=0.2)=h/8;
z0m(NDVI<0.2)=0.0005;
z0m=z0m';
% L=(zm-d0)./Lz;
Tf=(Ta+T0).*(uf.^2)./(k*g*L);
z0h=(70*nu./uf).*exp(-7.2*((uf.^0.5).*(abs(Tf).^0.25)));

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
         Fm(i,1)=NaN;
         Fh(i,1)=NaN;
        x(i,1)=NaN;
        x0(i,1)=NaN;
        y(i,1)=NaN;
        y0(i,1)=NaN;
      end
end
Ra=log(zh./z0h-Fh).*log(zm./z0m-Fm)./(U*k*k);
Rss=exp(7.99-3.4*(swc/0.65));

LEs=(deta.*As+Cp*rho_a.*VPD./Ra).*(1-Fw)./(deta+gama.*(1+Rss./Ra));
for i=1:m
   LEc(i,1)=LE(i,1)-LEs(i,1);
end 
Rsc=Cp*rho_a.*VPD./(gama.*LEc)+((deta.*(Ac-LEc)./(gama.*LEc))-1)./Ra;

Gs=1.0./Rsc;

for i=1:m
       if Ta(i,1)>=T_open
           m_Ta(i,1)=1.0;
       end
       if (Ta(i,1)<T_open)&&(Ta(i,1)>T_close)
           m_Ta(i,1)=(Ta(i,1)-T_close)/(T_open-T_close);
       end
       if Ta(i,1)<=T_close
           m_Ta(i,1)=0.1;
       end
end

for i=1:m  
    if (VPD(i,1)<=VPD_open)
       m_VPD(i,1)=1.0;
    end
    if (VPD(i,1)<VPD_close)&&(VPD(i,1)>VPD_open)
       m_VPD(i,1)=(VPD_close-VPD(i,1))/(VPD_close-VPD_open); 
    end
    if (VPD(i,1)>=VPD_close)
       m_VPD(i,1)=0.1;
    end
end

CL=Gs./(m_Ta.*m_VPD.*LAI);
CM=mean(CL);
RMSE=0;
MB=0;
R1=0;
R2=0;
R3=0;
a=0;
z=0;

for k=0.0001:0.0001:0.015
        Gst=k*m_Ta.*m_VPD.*LAI;
%       Rsct=1.0./Gst;
        LEct=(deta.*Ac+Cp*rho_a.*VPD./Ra).*(1-Fw)./(deta+gama.*(1+1.0./(Ra.*Gst)));
        LEt=LEs+LEct;
        M1=mean(LE);
        M2=mean(LEt);
        a=a+1;
        for i=1:m
           MB=MB+(LE(i)-LEt(i));
        end
        C=corrcoef(LE,LEt);
        p=polyfit(LE,LEt,1);
        R(a,1)=k;
        R(a,2)=sqrt(mean((LE-LEt).^2));
        R(a,3)=MB/m;
        R(a,4)=C(1,2);
        R(a,5)=p(1);
        R(a,6)=p(2);
        R(a,7)=M1;
        R(a,8)=M2;
        RMSE=0;
        MB=0;
        R1=0;
        R2=0;
        GS(:,a)=Gst;
        LC(:,a)=LEct;
        LS(:,a)=LEt;
end
% 
% Fs=(swd/1000).*(1100+b1)./(swd+b1);

% FTa=1-0.0016*(25-Ta).^2;
% 
% FD=1-0.0108*VPD;
% 
% for i=1:m
%     if swc(i,1)>0.75*0.65
%          Fswc(i,1)=1;
%     elseif swc(i,1)<0.08;
%          Fswc(i,1)=0;
%     elseif swc(i,1)>0.08&&swc(i,1)<0.75*0.65
%        Fswc(i,1)=(swc(i,1)-0.75*0.65)/(0.08-0.75*0.65);
%     end
% end
% for i=1:m
%     if LAI(i,1)<=2
%         LAIa(i,1)=LAI(i,1);
%     elseif LAI(i,1)>=4
%         LAIa(i,1)=0.5*LAI(i,1);
%     elseif LAI(i,1)<=4&&LAI(i,1)>=2
%        LAIa(i,1)=2;
%     end
% end
% 
% Fk=Gs*150./(LAIa.*Fswc.*FD.*FTa);




