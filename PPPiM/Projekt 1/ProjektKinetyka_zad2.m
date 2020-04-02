clear
clc
close all

g=9.811;
ro=1100;
mip=1.1;
tau0=25;
L=4;
delta=7e-3;
p1=1.24e5;
p2=1.1e5;
alpha=20.*pi/180;
% g=9.811;
% ro=1200;
% mip=1.5;
% tau0=23;
% L=7;
% delta=6e-3;
% p1=1.35e5;
% p2=1.4e5;
% alpha=35.*pi/180;

y=0:5e-4:delta;
tauxy=((p1-p2)./L+ro.*sin(alpha).*g).*y;
dzban=interp1(tauxy,y,tau0);
tab=[y'*1e3 tauxy'];
save('naprezenia.txt','tab','-ascii','-double');

figure(1)
grid on;
grid minor
xlabel('Prêdkoœæ cieczy [m/s]');
ylabel('Odleg³oœæ od œrodka pomiêdzy p³ytami [m]');
title('Rozk³ad naprê¿eñ');
hold all

plot(tauxy,y,'r');
plot(tauxy,-y,'r');
saveas(gcf,'F:\Pulpit\Nauka\Studia\III rok\PPPiM\Projekt1\TemperatureMeasurement_SCS\bin_nap','epsc')


figure(2)
grid on;
grid minor
xlabel('Prêdkoœæ cieczy [m/s]');
ylabel('Odleg³oœæ od œrodka pomiêdzy p³ytami [m]');
title('Rozk³ad prêdkoœci dla p³ynu binghamowskiego');
hold all

y=3.5e-3:2.5e-4:delta;
vx=(((p1-p2)./L+ro.*sin(alpha).*g).*(y.^2-delta.^2)./2-tau0.*(y-delta))./mip;
dzban2=min(vx);
tab2=[y'*1e3 -vx'];
save('predkosc.txt','tab2','-ascii','-tabs');

plot(-vx,y,'r');
plot(-vx,-y,'r');
plot(-[dzban2 dzban2],[dzban -dzban],'r');
saveas(gcf,'F:\Pulpit\Nauka\Studia\III rok\PPPiM\Projekt1\TemperatureMeasurement_SCS\bin_v','epsc')

Qst=-dzban2*dzban;
Qcal=integral(@(y) -(((p1-p2)./L+ro.*sin(alpha).*g).*(y.^2-delta.^2)./2-tau0.*(y-delta))./mip,dzban,delta);
Qcalanal=(((p1-p2)./L+ro.*sin(alpha).*g).*(dzban.^3-delta.^2)./6-tau0./2.*(dzban.^2-delta))./mip-...
    (((p1-p2)./L+ro.*sin(alpha).*g).*(delta.^3-delta.^2)./6-tau0./2.*(delta.^2-delta))./mip;
Qsum=Qst+Qcalanal;
m=Qsum*ro;














