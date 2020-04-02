%% Pawe³ Antkowiak
%Projekt 4 
%Aparatura

clear
clc
close all
%% Obliczenia

dw=1;%[m]
dn=0;%[m]
p0=8e5;%[Pa]
kr=240e6;%[Pa]
c1=0.002;
c2=0.001;
c3=0.0008;
hw=0.4;%[m]


err=1;
g(1)=8.5e-3;
i=1;
while err>1e-5
dz(i)=dw+2*g(i);
hz(i)=hw+g(i);
yw(i)=wspz(hz(i));
g1(i)=p0*dz(i)*yw(i)/4/kr;
err=abs(g(i)-g1(i));
g(i+1)=g1(i);
i=i+1;
end
function yw=wspz(om)
om1=[0.18 0.2 0.25 0.3 0.35 0.4 0.45 0.5];
yw1=[3.37 2.9 2 1.53 1.3 1.18 1.12 1.1];
yw=interp1(om1,yw1,om,'PCHIP');
end