% In¿ynieria reaktorów chemicznych, projekt 3, zadanie 2.
clc
clear
%% Dane
% Wspólne - wszystko oprócz ciœnienia w jednostkach podstawowych SI
kgA = 10;
klA = 3e-5;
HA = 125e-6;
L = 200;
G = 30;
C0 = 55555; %[mol/m3]

% Wiktor
pa2 = 1.1e-2;   %[bar]
pa1 = 9e-5;     %[bar]
DaDb = 0.8;

% Kasia
% pa2 = 1.5e-2;   %[bar]
% pa1 = 8e-5;     %[bar]
% DaDb = 0.6;

% Eliza
% pa2 = 1e-2;   %[bar]
% pa1 = 3e-5;   %[bar]
% DaDb = 0.6;

%% Obliczenia
% 1) Obliczyæ cBkryt
cBkryt = kgA*DaDb*pa2/klA

%% 2)
% cB1 - stê¿enie na wlocie do wie¿y
cB1 = cBkryt-((G*C0/L)*(pa1-pa2))

% Wysokoœæ kalumni
KGa = 1/((1/kgA)+(HA/klA));
h1 = (G/(KGa*(1-HA*G*C0/(DaDb*L))))*log((pa2+((pa1-pa2)*HA*G*C0/(DaDb*L))+HA*cBkryt/DaDb)/(pa1+HA*cBkryt/DaDb))

%% 3)
% Przep³yw cieczy
L = G*C0*(pa2-pa1)/cBkryt

% Wysokoœæ kalumni
h2 = (G/(KGa*(1+HA*G*C0/(DaDb*L))))*log((pa2+(pa2-pa1)*HA*G*C0/(DaDb*L))/pa1)