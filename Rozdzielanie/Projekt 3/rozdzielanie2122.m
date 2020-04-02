%% Pawe³ Antkowiak
clear
close all
clc

%% dane
xs=0.41;
xd=0.97;
xw=0.013;
alfa=2.43;
a=1.25;


%% Obliczenia
x1=0:0.1:1;
y1=x1;
y2=alfa.*x1./(1+(alfa-1).*x1);

xw1=0:1e-2:1;
yw1=xw1;
yw2=interp1(x1,y2,xw1,'PCHIP');

% plot1(xw1,yw1,1,'r')
% plot(xw1,yw2,'r')

ys=interp1(x1,y2,xs,'PCHIP');
% plot(xs,ys,'bo','LineWidth',2)
% plot(xd,xd,'bo','LineWidth',2)
% plot([xs xs],[0 1],'b')
% plot([xd xs],[xd ys],'b')
% text(0.55,0.7,'GLO');
% text(xs,0.9,'LP');

a2=(ys-xd)/(xs-xd);
b2=ys-a2*xs;

rmin=a2/(1-a2);

plot1(xw1,yw1,2,'r')
plot(xw1,yw2,'r')
plot([xd xd],[0 xd],'b');

schodki(xw1,yw2,xd,xw)
r=rmin*a;

a3=r/(r+1);
b3=1/(r+1);

ys2=a3.*xs+b3*xd;
plot1(xw1,yw1,3,'r')
plot(xw1,yw2,'r')
plot([xs xs],[0 1],'b')
plot(xs,ys2,'bo','LineWidth',2)
plot(xd,xd,'bo','LineWidth',2)
plot(xw,xw,'bo','LineWidth',2)
plot([xd xs],[xd ys2],'b')
plot([xw xs],[xw ys2],'b')

a4=(ys2-xw)/(xs-xw);
b4=ys2-a4*xs;

nt=schodki2(xw1,yw2,xd,xw,xs,a3,b3,a4,b4);

%% anale
rminu=1./(alfa-1).*(xd/xs-alfa.*(1-xd)./(1-xs));
nminf=log10(xd.*(1-xw)./(1-xd)./xw)./log10(alfa);

m=a3;
b5=1./(r+1).*xd;
k=fsolve(@(k) m.*(alfa-1).*k.^2+(m+b5.*(alfa-1)-alfa).*k+b5,0.5);
c=1+(alfa-1).*k;
M=m.*c.*(alfa-1)./(alfa-m.*c.^2);

gora=(xd-k).*(1-M.*(xs-k));
dol=(xs-k).*(1-M.*(xd-k));
dold=alfa./m./c.^2;
ng=log10(gora./dol)./log10(dold);

md=(ys-xw)/(xs-xw);
b5d=xw.*(1-md);
kd=fsolve(@(k) md.*(alfa-1).*k.^2+(md+b5d.*(alfa-1)-alfa).*k+b5d,0);
cd=1+(alfa-1).*kd;
Md=md.*cd.*(alfa-1)./(alfa-md.*cd.^2);

gora=(xs-kd).*(1-Md.*(xs-kd));
dol=(xw-kd).*(1-Md.*(xd-kd));
dold=alfa./md./cd.^2;
nd=log10(gora./dol)./log10(dold)

%% Funkcje
function schodki(x,y,xd,xw)

dzban1=1;
dzban2(1)=xd;
dzban3(1)=xd;
i=2;
% dzban2(2)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban2(1),dzban2(1));
% dzban3(2)=interp1(x,y,dzban2(2),'PCHIP');
% dzban2(3)=dzban2(2);
% dzban3(3)=dzban2(2);
while dzban1>xw
    dzban2(i)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban2(i-1),dzban2(i-1));
    dzban3(i)=interp1(x,y,dzban2(i),'PCHIP');
    dzban2(i+1)=dzban2(i);
    dzban3(i+1)=dzban2(i);
    dzban1=dzban2(i+1);
    i=i+2;
end
plot(dzban2,dzban3,'b');
end

function odp=schodki2(x,y,xd,xw,xs,a3,b3,a4,b4)

dzban1=1;
dzban2(1)=xd;
dzban3(1)=xd;
i=2;
% dzban2(2)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban2(1),dzban2(1));
% dzban3(2)=interp1(x,y,dzban2(2),'PCHIP');
% dzban2(3)=dzban2(2);
% dzban3(3)=dzban2(2);
while dzban1>xs
    dzban2(i)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban3(i-1),dzban2(i-1));
    dzban3(i)=interp1(x,y,dzban2(i),'PCHIP');
    dzban2(i+1)=dzban2(i);
    dzban3(i+1)=a3.*dzban2(i)+b3*xd;
    dzban1=dzban2(i+1);
    i=i+2;
end
i=i-1;
while dzban1>xw
    dzban2(i)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban3(i-1),dzban2(i-1));
    dzban3(i)=interp1(x,y,dzban2(i),'PCHIP');
    dzban2(i+1)=dzban2(i);
    dzban3(i+1)=a4.*dzban2(i+1)+b4;
    dzban1=dzban2(i+1);
    i=i+2;
end
plot(dzban2,dzban3,'b');
odp=max(size(dzban2))/2;
end














