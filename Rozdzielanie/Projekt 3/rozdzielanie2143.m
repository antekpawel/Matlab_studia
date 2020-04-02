%% Pawe³ Antkowiak
clear
close all
clc

%% dane
p=1.*101300;
xz=0.32;xz=xz./46./(xz./46+(1-xz)./18);
xs=0.54;xs=xs./46./(xs./46+(1-xs)./18);
xl=0.65;xl=xl./46./(xl./46+(1-xl)./18);
xd=0.84;xd=xd./46./(xd./46+(1-xd)./18);
xw=0.12;xw=xw./46./(xw./46+(1-xw)./18);
sd=1.80;
ld=0.46;
r=2.2;
rd=1.156;

%% Obliczenia

F=fsolve(@(x) dzban(xz,xs,xl,xd,xw,sd,ld,x),[1 1]);
wd=F(2);
zd=F(1);

i1=(r+1).*rd;
x2=(xd+ld*xl)/(1+ld);
x3=(xd+ld*xl-sd*xs)/(1+ld-sd);

%% Funkcje
function f=dzban(xz,xs,xl,xd,xw,sd,ld,x)
f(1)=sd+x(1)-1-x(2)-ld;
f(2)=sd.*xs+x(1).*xz-xd-x(2).*xw-ld.*xl;
end




