%% Pawe³ Antkowiak
clear
close all
clc

%% dane
t=20+273;
P=2.9e5;
m=0.15;
x0=0.01;
y1=120e-4;
y0=40e-4;
a=1.44;
N=9;
n=1.1;

%% Obliczenia
%a
ymax=0.1;
x1max=y1./m;
ldogmin=(y0-y1)./(x0-x1max);
plot1([0 ymax],[0 m.*ymax],1,'r')
plot([x0 x1max],[y0 x1max.*ldogmin+y0-ldogmin*x0],'b')
plot([0 x1max x1max],[y1 y1 0],'k--');
%b
ldog=a*ldogmin;
x1=(y1+ldog*x0-y0)/ldog;
trojzabposejdona=(y1-y0)/(y1-m*x0);
A=ldog/m;
Nob=log10((trojzabposejdona-A)/(trojzabposejdona-1))/log10(A)-1;
%c
mp=y1/x1;
xsr=(x0+x1)/2;%pomaga przy H
% dowiedziec sie co z H
H=1.8e3/xsr;
p0=H/mp*18/29;
%d
yp0=y1-(A^(N+1)-A)/(A^(N+1)-1)*(y1-m*x0);
%e
xp1max=yp0/m;
nmax=1/ldog*(y1-yp0)/(xp1max-x0);
ldogr=n*ldog;
x0max=x1-ldogr*(y1-y0);
Ap=ldogr/m;
trojzabposejdonap=2.5/((y1-yp0)/(y1-m*x0max));
Nobp=real(log10((trojzabposejdonap-Ap)/(trojzabposejdonap-1))/log10(Ap)-1);




