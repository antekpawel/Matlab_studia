%% Pawe³ Antkowiak
clear
close all
clc

%% dane

ya1=0.033;
yb1=0.058;
yc1=0.040;
t=25+273;
alfa=1.4;
eta=0.95;

%odczyane
ka=1396;
kb=338;
kc=1012;

%% Obliczenia

yb0=(1-eta)*yb1;
xb1=yb1/kb;
ldogmin=(yb1-yb0)/xb1;
ldog=alfa*ldogmin;
psib=(yb1-yb0)/yb1;
ab=ldog/kb;
N=log10((psib-ab)/(psib-1))/log10(ab)-1;
ya0=ya1-((ldog/ka)^(N+1)-(ldog/ka))/((ldog/ka)^(N+1)-1)*ya1;
yc0=yc1-((ldog/kc)^(N+1)-(ldog/kc))/((ldog/kc)^(N+1)-1)*yc1;
etaa=(ya1-ya0)/ya1
etac=(yc1-yc0)/yc1



