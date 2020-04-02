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
kr=240e6/1.65;%[Pa]
c1=0.002;
c2=0.001;
c3=0.0008;


err=1;
g(1)=15e-3;
i=1;
while err>1e-5
dz(i)=dw+2*g(i);
beta(i)=dz(i)/dw;
om(i)=dn/sqrt(dz(i)*g(i));
z(i)=wspz(om(i));
g1(i)=p0*dz(i)/(2.3*kr*z(i)+p0);
err=abs(g(i)-g1(i));
g(i+1)=g1(i);
i=i+1;
end

function z=wspz(om)
om1=[0 0.5 1 1.5 2 3 4 5];
z1=[1 0.8 0.64 0.53 0.44 0.32 0.24 0.18];
z=interp1(om1,z1,om,'PCHIP');
end









