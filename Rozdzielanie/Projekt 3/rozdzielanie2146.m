%% Pawe³ Antkowiak
clear
close all
clc

%% Dane
n=3;
aab=3;
r=2;
k=2;
x0=0.5;
a=0.6;
% n=5;
% aab=2.8;
% r=2.0;
% k=2;
% x0=0.5;
% a=0.38;
%% Obliczenia

%zakladam xd

xd=0.9003;
xw=0:1e-2:1;
yw=xw;

plot1(xw,yw,1,'r');
yw2=aab.*xw./(1+(aab-1).*xw);
plot(xw,yw2,'r')
a1=r/(r+1);
b1=1/(r+1).*xd;
plot([0.1 xd],[0.1.*a1+b1 xd],'b');
xm=schodki1(xw,yw2,xd,n,k,a1,b1);
ym=xm.*a1+b1;
xdm=(xd+xm)/2;
plot([xdm xdm],[0 xdm],'g')
plot([xd xd],[0 xd],'g')
plot([xm xm],[0 ym],'g')

a2=(xdm-ym)/(xdm-xm);
b2=ym-a2*xm;
plot([0.2 xd],[0.2*a2+b2 xd*a2+b2],'b')
x1=schodki2(xw,yw2,xm,n,k,a2,b2);
y1=a2*x1+b2;
plot([x1 x1],[0 y1],'g')

%% tera jebaj z tego tablice


wek=1:1:15;
for i=wek
xd(i)=0.01.*i+0.75;

b1(i)=1/(r+1).*xd(i);
xm(i)=schodki3(xw,yw2,xd(i),n,k,a1,b1(i));
ym(i)=xm(i).*a1+b1(i);
xdm(i)=(xd(i)+xm(i))/2;
a2(i)=(xdm(i)-ym(i))/(xdm(i)-xm(i));
b2(i)=ym(i)-a2(i)*xm(i);
x1(i)=schodki4(xw,yw2,xm(i),n,k,a2(i),b2(i));
y1(i)=a2(i)*x1(i)+b2(i);
end
tab1(:,1)=xd;
tab1(:,2)=xm;
tab1(:,3)=xdm;
tab1(:,4)=x1;

clc
for i=wek
ldol(i)=exp(integral(@(g) (xdm(i)-g).^-1,x0,x1(i)));
end
tab1(:,5)=ldol;
xdk=fsolve(@(g) interp1(xd,ldol,g)-1+a,0.9);
xk=interp1(xd,x1,xdk);
xmk=interp1(xd,xm,xdk);

dl=(1-ldol)/2;
tab1(:,6)=dl;
dk=interp1(xd,dl,xdk);

xdsr=(xdk+x0)./2
xmsr=(xmk+x0)./2

%% Funkcje 
% function f=fdzban(xd,a,a1,r,xw,yw2,n,k,x0)
% b1=1/(r+1).*xd;
% 
% [xdm,xm]=xdm1(r,xd,xw,yw2,n,k,a1)
% ym=xm.*a1+b1;
% a2=(xdm-ym)/(xdm-xm);
% b2=ym-a2*xm;
% plot([0.2 xd],[0.2*a2+b2 xd*a2+b2],'b')
% x1=schodki2(xw,yw2,xm,n,k,a2,b2);
% 
% pom=integral(@(x) 1./(xdm-x),x0,x1);
% f=1-a-exp(pom);
% 
% 
% 
% end
% function [xdm,xm]=xdm1(r,xd,xw,yw2,n,k,a1)
% b1=1/(r+1).*xd;
% xm=schodki3(xw,yw2,xd,n,k,a1,b1);
% xdm=(xd+xm)/2;
% end




function odp=schodki1(x,y,xd,n,k,a1,b1)

dzban2(1)=xd;
dzban3(1)=xd;

for i=2:2:2*(n-k)
    dzban2(i)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban3(i-1),dzban2(i-1));
    dzban3(i)=interp1(x,y,dzban2(i),'PCHIP');
    dzban2(i+1)=dzban2(i);
    dzban3(i+1)=a1.*dzban2(i)+b1;
    dzban1=dzban2(i+1);
end

plot(dzban2,dzban3,'b');
odp=min(dzban2);
end

function odp=schodki2(x,y,xd,n,k,a1,b1)

dzban2(1)=xd;
dzban3(1)=xd.*a1+b1;

for i=2:2:2*k
    dzban2(i)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban3(i-1),dzban2(i-1));
    dzban3(i)=interp1(x,y,dzban2(i),'PCHIP');
    dzban2(i+1)=dzban2(i);
    dzban3(i+1)=a1.*dzban2(i)+b1;
    dzban1=dzban2(i+1);
end

plot(dzban2,dzban3,'b');
odp=min(dzban2);
end
function odp=schodki3(x,y,xd,n,k,a1,b1)

dzban2(1)=xd;
dzban3(1)=xd;

for i=2:2:2*(n-k)
    dzban2(i)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban3(i-1),dzban2(i-1));
    dzban3(i)=interp1(x,y,dzban2(i),'PCHIP');
    dzban2(i+1)=dzban2(i);
    dzban3(i+1)=a1.*dzban2(i)+b1;
    dzban1=dzban2(i+1);
end

odp=min(dzban2);
end

function odp=schodki4(x,y,xd,n,k,a1,b1)

dzban2(1)=xd;
dzban3(1)=xd.*a1+b1;

for i=2:2:2*k
    dzban2(i)=fsolve(@(d) interp1(x,y,d,'PCHIP')-dzban3(i-1),dzban2(i-1));
    dzban3(i)=interp1(x,y,dzban2(i),'PCHIP');
    dzban2(i+1)=dzban2(i);
    dzban3(i+1)=a1.*dzban2(i)+b1;
    dzban1=dzban2(i+1);
end

odp=min(dzban2);
end
