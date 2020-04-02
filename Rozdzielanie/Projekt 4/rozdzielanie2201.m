%% Pawe³ Antkowiak
clear
close all
clc

%% dane

y1=29e-3;
eta=0.92;
x0=0.8e-3;
g=530/3600;
kga=4.1e-3*1e-5;
kca=3.05e-4;
a1=2.35;
a2=1.9;
P=1e5;

% z tabeli 3.11 funkcja interpnakrzywyryj
% z p1
x1=0.058;
% z p2
x0p=0.0056;
%% Obliczenia

y0=(1-eta)*y1;
p1=P*y1/(1-y1)
p2=P*y0/(1-y0)
ldogmin=(y1-y0)/(x1-x0);
ldog=a1*ldogmin;
x1p=(y1-y0)./ldog+x0;

ldogmax=(y1-y0)/(x1p-x0p);
nmax=ldogmax/ldog;
ldogn=ldog*(1+(nmax-1)/a2);
a=P*29/17;
m=y1/x1;
h=58.8/a/m;

Kga=(1/kga+1/h/kca)^-1;
A=1-1/ldog*m;
B=m*(1/ldog*y0-x0);
htu=g/17/a/Kga/a2;
nog=1/A*log((A*y1+B)/(A*y0+B));
wys=htu*nog;








