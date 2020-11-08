Data=xlsread('E:\Pre_Data\Ali\Test\Ali3.xlsx');
[m,n]=size(Data);

beta=50;
T0=273.15;
P0=1013.25; 

rd=287.0;
h=0.03;
k=0.41;
g=9.8;
ksac=0.622;

zh=Data(:,4);
zm=zh;
d0=Data(:,5);
U=Data(:,6);
Ta=Data(:,7);
RH=Data(:,8);
P=Data(:,9)/10;
NDVI=Data(:,10);
Rn=Data(:,11);
G0=Data(:,12);
uf=Data(:,13);
L=Data(:,14);
LE=Data(:,15);

Cp=4.2*0.242;
lamada=2501-2.361*Ta; 
rho_a=3.846.*P./(Ta+273.15);
gama=Cp*P./(0.622*lamada);
Es=0.6108*exp((17.27*Ta)./(237.3+Ta));
Ea=Es.*RH/100;
VPD=Es-Ea;
deta=4098*Es./((Ta+237.3).^2);

fc=(NDVI-0.2)/(0.8-0.2);

nu=1.327e-5*(P0./P).*(((Ta+T0)/ T0).^1.754); 

for i=1:m;
   if NDVI(i,1)<0.3
      Ac(i,1)=0;
      As(i,1)=Rn(i,1)-G0(i,1);
      LEs(i,1)=LE(i,1);
   else
      Ac(i,1)=Rn(i,1)*fc(i,1);
      As(i,1)=Rn(i,1)*(1-fc(i,1))-G0(i,1);
      LEs(i,1)=NaN;
   end
end


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

z0m=0.138*h;

%L=(zm-d0)./Lz;
Tf=Ta.*(uf.^2)./(k*g*L);
z0h=(70*nu./uf).*exp(-7.2*((uf.^0.5).*(abs(Tf).^0.25)));

for i=1:m
  if isnan(L(i,1))==0
     if L(i,1)>0
        Fm(i,1)=-5.3*(zm(i,1)-z0m)/L(i,1);
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
          x0(i,1)=(1-19*z0m/L(i,1))^0.25;
          y(i,1)=(1-11.6*zh(i,1)/L(i,1))^0.5;
          y0(i,1)=(1-11.6*z0h(i,1)/L(i,1))^0.5;
          Fm(i,1)=2*log((1+x(i,1))/(1+x0(i,1)))+log((1+x(i,1)*x(i,1))/(1+x0(i,1)*x0(i,1)))-2*atan(x(i,1))+2*atan(x0(i,1));
          Fh(i,1)=2*log((1+y(i,1))/(1+y0(i,1)));
       end
end
end
Ra=log(zh./z0h-Fh).*log(zm/z0m-Fm)./(U*k*k);

for i=1:m
    if NDVI(i,1)<0.3  
         LEs(i,1)=LE(i,1);
    else
         LEs(i,1)=NaN;
    end
end

Rs=(((deta.*As+Cp*rho_a.*VPD./Ra).*(1-Fw).*Fs./LE-deta)./gama-1).*Ra;










