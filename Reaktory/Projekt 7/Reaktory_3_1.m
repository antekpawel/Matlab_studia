%% Pawe³ Antkowiak
clear
close all
clc

%% Dane

car0=1300;
qr=6e-4;
fi=0.07;
ma=1.6;
k1=5;
k2=0.25;
cbr0=0;
cae0=0;
cbe0=0;
alfaa=0.95;

%% Obliczenia

qe=qr/fi; 
car=fi*car0*(1-alfaa)/(fi+ma);
cae=car*ma;

carp=car*qr/(qr+qe);
caep=cae*qe/(qr+qe);
car0p=car0*qr/(qr+qe);

qc=qe+qr;
cbep=car0p-caep-carp;
tau=cbep/(k1*caep-k2*cbep);

v=qc*tau;



