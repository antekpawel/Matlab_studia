%% Pawe³ Antkowiak
%Projekt 2 zad2
%Kinetyka

clear
clc
close all
K=273.15;
g=9.81;
%% Dane projektowe

d=0.1e-3;
tw=40+K;
tp=55+K;
p=101315;
D=2.636e-5;
dk=0.02e-3;

%% Dane odczytane
tsr=(tw+tp)./2;

row=989.20;
miw=5.719e-4;
niw=miw./row;
Mw=18.015e-3;

rop=1.1010;
mip=1.9487e-5;
nip=1.7695e-5;
Mp=29e-3;

pnas=39566;
R=8.314;
g=9.81;

Dsr=D.*(tsr./(25+K)).^1.5;
beta=(Mp-Mw)./rop;




%% Obliczenia prêdkoœci
u1=d.^2.*(row-rop).*g./18./mip;
u2=0.153.*d.^1.14.*(row-rop).^0.71.*g.^0.71./rop.^0.29./mip.^0.43;
u3=1.74.*(d.*(row-rop).*g./rop).^0.5;

Re1=u1.*rop.*d./mip;
Re2=u2.*rop.*d./mip;
Re3=u3.*rop.*d./mip;

Gr=g.*beta.*d.^3.*pnas./R./tsr./niw.^2;
Sc=mip./rop./Dsr;
disp(Gr.*Sc);
Sh1=2+0.282.*(Gr.*Sc).^0.37;%GrSc<100
Sh2=2+0.5.*(Gr.*Sc).^0.25;
Sh3=2+0.6.*Re2.^0.5.*Sc.^0.33;
kp=Sh3.*Dsr./d./R./tsr;
Na=kp.*(pnas).*Mw.*4.*pi.*(d./2).^2;
t=integral(@(x) time(x),d,dk);
kp2=Sh1.*Dsr./d./R./tsr;
w=integral(@(x) speed(x),dk,d);
t
w
m1=row*pi*d^3/6;
m2=row*pi*dk^3/6;
t1=(m1-m2)/w

a=2*Dsr*pnas*Mw/R/tsr/row;
b=(0.1*Sc^0.33*sqrt((row-rop)*rop*g/2/mip^2))^(-0.333);
x=sqrt(d/2);
t2=b^4/a*(x/b-log((1+x/b)^2/(1-x/b+(x/b)^2))/6-(atan((2*x/b-1)/sqrt(3))+atan(1/sqrt(3)))/sqrt(3))
x=sqrt(dk/2);
t3=b^4/a*(x/b-log((1+x/b)^2/(1-x/b+(x/b)^2))/6-(atan((2*x/b-1)/sqrt(3))+atan(1/sqrt(3)))/sqrt(3))
t=t2-t3
m0=row/6*pi*d^3
mk=row/6*pi*dk^3
chuj=(m0-mk)/w
h=integral(@(x) high(x),dk,d)









%% Funkcje u¿ywane w rozwi¹zaniu problemu

%Funkcja do obliczania gêstoœci wody na podstawie danych z NIST'a
function ro=Ro(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
RO=I(:,3);
ro=interp1(T,RO,t,'PCHIP');
end

%Funkcja do obliczania lepkoœci dynamicznej wody na podstawie danych z NIST'a
function mi=Mi(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
RO=I(:,12);
mi=interp1(T,RO,t,'PCHIP');
end

function t=time(d)
K=273.15;
tw=40+K;
tp=55+K;
D=2.636e-5;
tsr=(tw+tp)./2;
row=989.20;
Mw=18.015e-3;
rop=1.1010;
mip=1.9487e-5;
pnas=39566;
R=8.314;
g=9.81;
Dsr=D.*(tsr./(25+K)).^1.5;
Sc=mip./rop./Dsr;

u2=0.153.*d.^1.14.*(row-rop).^0.71.*g.^0.71./rop.^0.29./mip.^0.43;
Re=d.*u2.*rop./mip;
Sh=2+0.6.*Re.^0.5.*Sc.^0.33;
licznik=-row.*R.*tsr.*d;
mianownik=Sh.*2.*D.*pnas.*Mw;
t=licznik./mianownik;
end

function w=speed(d)
K=273.15;
rop=1.1010;
mip=1.9487e-5;
nip=1.7695e-5;
row=989.20;
miw=5.719e-4;
tw=40+K;
tp=55+K;
D=2.636e-5;
tsr=(tw+tp)./2;
Mw=18.015e-3;
rop=1.1010;
mip=1.9487e-5;
niw=miw./row;
Mp=29e-3;
pnas=39566;
R=8.314;
g=9.81;
Dsr=D.*(tsr./(25+K)).^1.5;
beta=(Mp-Mw)./rop;
Gr=g.*beta.*d.^3.*pnas./R./tsr./niw.^2;
Sc=mip./rop./Dsr;
disp(Gr.*Sc);
Sh1=2+0.282.*(Gr.*Sc).^0.37;%GrSc<100


w=Dsr*pi./R./tsr.*Sh1.*pnas.*Mw;
end

function t=high(d)
K=273.15;
tw=40+K;
tp=55+K;
D=2.636e-5;
tsr=(tw+tp)./2;
row=989.20;
Mw=18.015e-3;
rop=1.1010;
mip=1.9487e-5;
pnas=39566;
R=8.314;
g=9.81;
Dsr=D.*(tsr./(25+K)).^1.5;
Sc=mip./rop./Dsr;

u2=0.153.*d.^1.14.*(row-rop).^0.71.*g.^0.71./rop.^0.29./mip.^0.43;
Re=d.*u2.*rop./mip;
Sh=2+0.6.*Re.^0.5.*Sc.^0.33;
licznik=u2;
mianownik=Sh.*2.*D.*pnas.*Mw./row./R./tsr./d;
t=licznik./mianownik;
end






