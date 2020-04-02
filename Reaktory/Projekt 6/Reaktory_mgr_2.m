%% Pawe³ Antkowiak
%Reaktory mgr. Projekt 1

clear
close all
clc

R=8.314;
K=273.15;

%% Dane projetowe

tau=100*60;
c0=1;

kg=2e-7*1e-4*1e3/60;
kn=4e4*1e-3*1e6/60;
cg=0.05;
rhok=4500;
ka=6;
kv=1;

N0=1;
v0=1;
dt1=tau;
a0=0.5;

%% Obliczenia zad1

c=c0;
err=1;

while abs(err)>1e-10
G=kg*(c-cg);
Rn=kn*(c-cg).^2;
n0=Rn/G;

m0=Rn*tau;
Nt=m0;
m1=n0*tau^2*G^2;
Lt=m1;
m2=2*m1*tau*G;
At=m2*ka;
m3=3*m2*tau*G;
Mt=m3*kv*rhok;
m4=4*m3*tau*G;
c1=c0-Mt;
err=c-c1;
c=c1;
end

L10=m1/m0;
L30=(m3/m0)^(1/3);
L21=m2/m1;
L32=m3/m2;
L43=m4/m3;

%% Obliczenia zad2

t0=0;
t=t0:1:t0+dt1;

teta=N0.*a0.*t;
M0=2./(2+teta);
M1=M0.^(2/3);
M2=M0.^(1/3);
M3=M0.^0;
M4=M0.^(-1/3);

figure(1)
clf

plot(t,M0,'r');
hold all
plot(t,M1,'b');
plot(t,M2,'g');
plot(t,M3,'k');
plot(t,M4,'m');

grid on
grid minor
axis([0 6000 0 12])

ts=fsolve(@(t) (2./(2+N0.*a0.*t)).^(-1/3)-2.*max(M4),6000)
ts/3600

t=logspace(0.1,log10(dt1),15);

teta=N0.*a0.*t;
M0=2./(2+teta);
M1=M0.^(2/3);
M2=M0.^(1/3);
M3=M0.^0;
M4=M0.^(-1/3);


plot(t,M0,'rx');
plot(t,M1,'bx');
plot(t,M2,'gx');
plot(t,M4,'mx');

tab1(:,1)=t;
tab1(:,2)=M0;
tab1(:,3)=M1;
tab1(:,4)=M2;
tab1(:,5)=M3;
tab1(:,6)=M4;

legend('M_0','M_1','M_2','M_3','M_4','Location','SouthWest')
podpis('Czas [s]','Moment znormalizowany')












