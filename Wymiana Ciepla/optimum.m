function mini=optimum
dz=0.14*1.075;
di=dz:1e-4:0.4;
k=koszt(di);
df(di);

mini=fsolve(@(di) df(di),0.25);

wykresy;
Ti=Tizo(mini)

function a=alfa(di,Ti)

T=[-50;-40;-30;-20;-10;0;10;20;30;40;50;60;70;80;90;100];
T=T+273.15;
l=[2.04000000000000;2.12000000000000;2.20000000000000;2.28000000000000;2.36000000000000;2.44000000000000;2.51000000000000;2.59000000000000;2.67000000000000;2.76000000000000;2.83000000000000;2.90000000000000;2.97000000000000;3.05000000000000;3.13000000000000;3.21000000000000];
l=l.*1e-2;
ni=[9.23000000000000;10.0400000000000;10.8000000000000;11.7100000000000;12.4300000000000;13.2800000000000;14.1500000000000;15.0600000000000;16;16.9600000000000;17.9500000000000;18.9700000000000;20.0200000000000;21.0900000000000;22.1000000000000;23.1300000000000];
ni=ni.*1e-6;

g=9.81665;
b=1./280;

a=0.5.*((g.*b.*(Ti-280).*spline(T,l,280).^4)./(di.*(spline(T,ni,280)).^2)).^(0.25);

function q=strumien(di)

Tw=346;
Tz=280;

Ti=Tizo(di);

dw=0.14;
dz=1.075.*dw;

li=0.009;
ls=50.3;

R1=(1./(2*pi*ls)).*log(dz./dw);
R2=(1./(2*pi*li)).*log(di./dz);
R3=1./(alfa(di,Ti).*pi.*di);
R=R1+R2+R3;

q=(Tw-Tz)./R;

function k=koszt_i(di)

dz=1.075.*0.14;
z=1507;
k=0.25.*pi.*(di.^2-dz.^2).*z;

function k=koszt_e(di)

zb=1.7e-8;
t=5*365*24*3600;
k=strumien(di).*zb.*t;

function k=koszt(di)

ke=koszt_e(di);
ki=koszt_i(di);
k=ke+ki;

function dff=df(di)

h=1e-6;
dff=(koszt(di-h)-koszt(di+h))./(2.*h);



function wykresy

T=[-50;-40;-30;-20;-10;0;10;20;30;40;50;60;70;80;90;100];
T=T+273.15;
l=[2.04000000000000;2.12000000000000;2.20000000000000;2.28000000000000;2.36000000000000;2.44000000000000;2.51000000000000;2.59000000000000;2.67000000000000;2.76000000000000;2.83000000000000;2.90000000000000;2.97000000000000;3.05000000000000;3.13000000000000;3.21000000000000];
l=l.*1e-2;
ni=[9.23000000000000;10.0400000000000;10.8000000000000;11.7100000000000;12.4300000000000;13.2800000000000;14.1500000000000;15.0600000000000;16;16.9600000000000;17.9500000000000;18.9700000000000;20.0200000000000;21.0900000000000;22.1000000000000;23.1300000000000];
ni=ni.*1e-6;

figure(1);

subplot(2,1,1);
plot(T,l,'b-',T,spline(T,l,T),'rx')
xlabel('T [K]')
ylabel('\lambda [W/mK]')
title('Interpolacja wspó³czynnika przewodzenia ciep³a w zakresie temperatur od -50^oC do 100^oC')
legend('Interpolacja metod¹ krzywych sklejanych','Dane literaturowe')
grid on

subplot(2,1,2);
plot(T,ni,'b-',T,spline(T,ni,T),'rx')
xlabel('T [K]')
ylabel('\nu [m^2/s]')
title('Interpolacja wspó³czynnika lepkoœci kinematycznej w zakresie temperatur od -50^oC do 100^oC')
legend('Interpolacja metod¹ krzywych sklejanych','Dane literaturowe')
grid on

dz=0.14*1.075;
di=dz:1e-4:0.4;

figure(2);

plot(di,koszt_e(di),'b-',di,koszt_i(di),'g-',di,koszt(di),'r-')
xlabel('Œrednica izolacji d_i [m]')
ylabel('Koszt [jedn. war]')
legend('Koszty eksploatacyjne','Koszty inwestycyjne','Koszt ca³kowity')
grid on
axis([dz 0.4 0 700]);


min=fsolve(@(di) df(di),0.25);

figure(3)

plot(di,koszt_e(di),'b-',di,koszt_i(di),'g-',di,koszt(di),'r-',min,koszt(min),'kx',[dz,min],[koszt(min),koszt(min)],'k--',[min,min],[koszt(min),0],'k--')
xlabel('Œrednica izolacji d_i [m]')
ylabel('Koszt [jedn. war]')
legend('Koszty eksploatacyjne','Koszty inwestycyjne','Koszt ca³kowity','Minimum kosztów')
grid on
axis([dz 0.4 0 700]);

function Ti=Tizo(di)

Tw=346;
Tp=280;
Ti=(Tp+Tw)./2;
t=0;
dw=0.14;
dz=1.075.*dw;
li=0.009;
ls=50.3;

while(abs(Ti-t)>1e-6)
    t=Ti;
    Q=(Tw-Ti)./(log(dz./dw)./(2.*pi.*ls)+log(di./dz)./(2.*pi.*li));
    Ti=Q./(pi.*di.*alfa(di,Ti))+Tp;
end

    


