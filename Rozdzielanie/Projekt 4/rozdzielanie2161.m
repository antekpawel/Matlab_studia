%% Pawe³ Antkowiak
clear
close all
clc

%% dane
xsa=0.47;
xda=0.94;
xwa=0.03;
a=1.16;

%% Obliczenia
rmin=fsolve(@(r) dzban(r,xsa,xda,xwa),5);

dz3=(xsa-xwa)./(xda-xwa);
xedmin=1./(dz3.*rmin+2);

r=rmin*a;
xed=1./(dz3.*r+2);
aabd=2.314.*xed+0.9214;

xeg=1./(dz3.*r+1);
aabg=2.314.*xeg+0.9214;
%% Rysowanie wykresu
x1=0:1e-2:1;
y1=x1;
x2=0:1e-3:xsa;
y2=x2.*aabd./(1+(aabd-1).*x2);
x3=xsa:1e-3:1;
y3=x3.*aabg./(1+(aabg-1).*x3);
plot1(x1,y1,1,'r');
plot(x2,y2,'r')
plot(x3,y3,'r')
plot([xsa xsa],[0.4 0.7],'b');

x4=x3;
y4=(r/(1+r)).*x4+(1/(1+r)).*xda;
plot(x4,y4,'b')

x5=[xwa xsa];
y5=[xwa (r/(1+r)).*xsa+(1/(1+r)).*xda];
plot(x5,y5,'b')

plot([xda xda],[0 xda],'b')

odp1=schodki2(x3,y3,xda,xsa,r/(1+r),(1/(1+r)).*xda);

ys2=(r/(1+r)).*xsa+(1/(1+r)).*xda;
a4=(ys2-xwa)/(xsa-xwa);
b4=ys2-a4*xsa;
odp2=schodki3(x2,y2,xsa,xwa,a4,b4);

%% Funkcje pomocnicze
function rmin=dzban(r,xsa,xda,xwa)

dz1=r./(r+1).*xsa+1./(r+1).*xda;
dz3=(xsa-xwa)./(xda-xwa);
xed=1./(dz3.*r+2);
aabd=2.314.*xed+0.9214;
dz2=aabd.*xsa./(1+(aabd-1).*xsa);
rmin=dz1-dz2;
end
%% Schodki
function odp=schodki2(x,y,xd,xs,a3,b3)

dzban1=1;
dzban2(1)=xd;
dzban3(1)=xd;
i=2;
% dzban2(2)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban2(1),dzban2(1));
% dzban3(2)=interp1(x,y,dzban2(2),'PCHIP');
% dzban2(3)=dzban2(2);
% dzban3(3)=dzban2(2);
while dzban1>=xs
    dzban2(i)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban3(i-1),dzban2(i-1));
    dzban3(i)=interp1(x,y,dzban2(i),'PCHIP');
    dzban2(i+1)=dzban2(i);
    dzban3(i+1)=a3.*dzban2(i)+b3;
    dzban1=dzban2(i+1);
    i=i+2;
end

plot(dzban2,dzban3,'b');
odp=max(size(dzban2))/2;
end



function odp=schodki3(x,y,xd,xs,a3,b3)

dzban1=1;
dzban2(1)=xd;
dzban3(1)=a3*xd+b3;
i=2;
% dzban2(2)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban2(1),dzban2(1));
% dzban3(2)=interp1(x,y,dzban2(2),'PCHIP');
% dzban2(3)=dzban2(2);
% dzban3(3)=dzban2(2);
while dzban1>=xs
    dzban2(i)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban3(i-1),dzban2(i-1));
    dzban3(i)=interp1(x,y,dzban2(i),'PCHIP');
    dzban2(i+1)=dzban2(i);
    dzban3(i+1)=a3.*dzban2(i)+b3;
    dzban1=dzban2(i+1);
    i=i+2;
end

plot(dzban2,dzban3,'b');
odp=max(size(dzban2))/2;
end









