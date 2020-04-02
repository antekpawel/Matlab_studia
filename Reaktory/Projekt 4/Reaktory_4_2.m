%% Nikt
clear
close all
clc

%% Dane

vr=10e-3;
k1=2e-4;
R=190e-6;
daeff=1e-10;
k1p=100*k1;
epsk=0.07;
ks=2.78e-7;

%% Obliczenia

%a
teta=R/3*sqrt(k1p./daeff);
ni=(3*teta*coth(3*teta)-1)/3./teta.^2;

q1=0:1e-6:1.5/3600;
alfa1=(k1.*(1-epsk).*vr+ni.*k1p.*epsk.*vr)./(q1+k1.*(1-epsk).*vr+ni.*k1p.*epsk.*vr);

wykres(1)
plot(q1*3600,alfa1,'b');
podpis('Przeplyw przez reaktor $$ [m^3/h]$$','Stopien przereagowania [-]');

q2=logspace(-5.5,-3.5,20);
alfa2=(k1.*(1-epsk).*vr+ni.*k1p.*epsk.*vr)./(q2+k1.*(1-epsk).*vr+ni.*k1p.*epsk.*vr);
tab(:,1)=q2*3600;
tab(:,2)=alfa2;

%b
k1pp=100*k1*exp(-ks*73*18e4);

q3=0:1e-6:1.5/3600;
alfa3=(k1.*(1-epsk).*vr+ni.*k1pp.*epsk.*vr)./(q3+k1.*(1-epsk).*vr+ni.*k1pp.*epsk.*vr);

plot(q3*3600,alfa3,'r');
podpis('Przeplyw przez reaktor $$ [m^3/h]$$','Stopien przereagowania [-]');

q4=logspace(-5.5,-3.5,20);
alfa4=(k1.*(1-epsk).*vr+ni.*k1pp.*epsk.*vr)./(q4+k1.*(1-epsk).*vr+ni.*k1pp.*epsk.*vr);
tab(:,3)=alfa4;


plot(q2*3600,alfa2,'bx');

plot(q4*3600,alfa4,'rx');

legend('Nowy katalizator','Katalizator po dwoch latach')

ile1=interp1(alfa1,q1,0.5)/interp1(alfa3,q3,0.5);

%c
teta2=R/3*sqrt(k1pp/daeff);
ni2=(3*teta2*coth(3*teta2)-1)/3./teta2.^2;










