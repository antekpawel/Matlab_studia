function k=koszt_e(di)

zb=1.7e-8;
t=5*365*24*3600;
k=strumien(di).*zb.*t;

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