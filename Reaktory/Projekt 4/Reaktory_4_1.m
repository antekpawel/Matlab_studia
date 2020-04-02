%% Pawe³ Antkowiak
clear
close all
clc

%% Dane

tzk=165;
R=190;
taup=130;

%% Obliczenia

A=taup/R^2;
rp=sqrt(tzk/A);
Rw=220:60:2500;
for i=1:max(size(Rw))
    rc(i)=fsolve(@(r) 2.*(1-r)-3.*(1-r).^(2/3)+1-tzk./(A.*Rw(i).^2),0.5);
end
clc

B=taup/R;
Rwb=fsolve(@(r) 1-(tzk./B./r).^3,200):30:2500;
alfab=1-(1-tzk./B./Rwb).^3;

wykres(1)
plot(Rw,rc,'r');
plot(Rwb,alfab,'b');

axis([rp 2500 0 1]);
podpis('Promien ziarna $$[\mu m]$$','Stopien przereagowania [-]');

Rt=300:100:2500;
for i=1:max(size(Rt))
    rct(i)=fsolve(@(r) -2.*r-3.*(1-r).^(2/3)+3-tzk./A./Rt(i).^2,0.5);
end
tab(:,1)=Rt;
tab(:,2)=rct;

Rtb=250:100:2500;
alfabtab=1-(1-tzk./B./Rtb).^3;
tab2(:,1)=Rtb;
tab2(:,2)=alfabtab;

plot(Rt,rct,'rx')
plot(Rtb,alfabtab,'bx')

legend('a) dyfuzja molekularna','b) reakcje powierzchniowa')      




