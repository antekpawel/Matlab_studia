%% Pawe³ Antkowiak & Magdalena Tomczak
%Projekt bioreaktora

clear
clc
close all
K=273.15;
g=9.81;

% xlswrite('C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Tabele\Tab1',proba);
% saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Wykresy\Wyk1','emf')
%% Dane projektowe

s0=17;%[kg/m^3]
mi_max=0.47;%[l/h]
ks=0.8;%[kg/m^3]
v=5.2;%[m^3]
yxs=0.55;%[-]

mi_max=mi_max/1000/3600;%przeliczenie na SI

%% Obliczenia

mi90=(mi_max/10*s0)/(ks+s0/10);
f90=mi90*v;

d90=f90/v;
x90=yxs*(s0-ks*d90/(mi_max-d90));
q90=d90*x90;

dopt=fsolve(@(d) yxs*(s0-(2*d*ks*(mi_max-d)+d^2*ks)/(mi_max-d)^2),0);
xopt=yxs*(s0-(ks*dopt)/(mi_max-dopt));
qopt=dopt*xopt;

porownanie=qopt/q90*100;





















