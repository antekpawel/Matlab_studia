%% Pawe³ Antkowiak
clear
close all
clc

%% Dane

kgaa=10;
klaa=3e-5;
ha=125e-6;
l=200;
g=30;
p0=1;
c0=55.555e3;

dadodb=0.8;
y2=1.1e-2;
y1=9e-5;

%% Obliczenia

%a
pa2=y2*p0;
pa1=y1*p0;
cbkr=pa2*kgaa/klaa*dadodb;

%b

h1=g*log(pa2/pa1)/p0/kgaa;
cb1=g*c0/p0/l*(pa2-pa1)+cbkr;

%c
cbkr2=pa1*kgaa/klaa*dadodb;
L=g*c0/p0/cbkr2*(pa2-pa1);
Kgaa = 1/((1/kgaa)+(ha/klaa));
h2=(g/(Kgaa*(1+ha*g*c0/(dadodb*L))))*log((pa2+(pa2-pa1)*ha*g*c0/(dadodb*L))/pa1);







