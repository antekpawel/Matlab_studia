%% PRCh Bilans Populacji - Zadanie 1
clc
clear

% Dane:
t = 100; %[min]
c0 = 3; %[g/dm^3]=[kg/m^3]
% Reszta przeliczona na SI
kg = (2e-7)/10;
kn = (4e+4)*1000;
cn = 0.05;
rok = 4.5e+3;
ka = 6;
kv = 1;

%% Obliczenia wstêpne
A = 6*rok*kv*(t^4)*(kg^3)*kn;
B = roots([A -5*A*cn 10*A*(cn^2) -10*A*(cn^3) (5*A*(cn^4))+1 (-A*(cn^5))-c0]);
c = B(5);
RN = kn*(c-cn)^2;
G = kg*(c-cn);
n0 = RN/G;

%% Momenty
m0 = RN*t;
m1 = n0*(t^2)*(G^2);
m2 = 2*n0*(t^3)*(G^3);
m3 = 6*n0*(t^4)*(G^4);
m4 = 24*n0*(t^5)*(G^5);

%% Wyniki
L10 = m1/m0;
L30 = (m3/m0)^(1/3);
L21 = m2/m1;
L32 = m3/m2;
L43 = m4/m3;
NT = m0 %[#/m^3]
AT = m2*ka %[m^2/m^3]
MT = m3*rok*kv %[kg/m^3]

%% Rozk³ad rozmiarów kryszta³ów
L = [0:0.1:35];
nL = n0*exp(-L./(G*t*(1e+6)));
figure(1)
hold all
grid minor
title('Rozk³ad rozmiarów kryszta³ów')
xlabel('L [{\mu}m]')
ylabel('n(L)')
plot(L,nL)