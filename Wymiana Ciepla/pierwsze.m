function teta=pierwsze

X=0.09;
%szk³o zwierciadlane
N=9;
Tk=359;
Tk=Tk+273.15;
Tp=396;
Tp=Tp+273.15;
Up=3;
Lcalk=22;
C=0.14;
n=0.75;

T0=22;
T0=T0+273.15;

teta=tbezw(Tp,Tk,T0);
a=alfa;

lszklo=0.884;

Bi1=Biot(a,lszklo,X./2);
Bi2=Biot(a,lszklo,N.*X./2);

f=fsolve(@(Fo) rownanie(Bi1,Bi2,N,X,Fo,0),800);
clc

l=lszklo;
x=0.09;
a=0.464.*10^6;
t=f.*(x.^2)./a


function teta=tbezw(Tp,Tk,T0)

teta=(Tk-Tp)./(T0-Tp);

function [l,ni]=wsp(T0,Tp)

t=(T0+Tp)/2;

T=[-50;-40;-30;-20;-10;0;10;20;30;40;50;60;70;80;90;100];
T=T+273.15;
l=[2.04000000000000;2.12000000000000;2.20000000000000;2.28000000000000;2.36000000000000;2.44000000000000;2.51000000000000;2.59000000000000;2.67000000000000;2.76000000000000;2.83000000000000;2.90000000000000;2.97000000000000;3.05000000000000;3.13000000000000;3.21000000000000];
l=l.*1e-2;
ni=[9.23000000000000;10.0400000000000;10.8000000000000;11.6100000000000;12.4300000000000;13.2800000000000;14.1500000000000;15.0600000000000;16;16.9600000000000;17.9500000000000;18.9700000000000;20.0200000000000;21.0900000000000;22.1000000000000;23.1300000000000];
ni=ni.*1e-6;

l=spline(T,l,t);
ni=spline(T,ni,t);

function a=alfa
[l,ni]=wsp(21+273.15,396+273.15);
c=0.14;
n=0.75;
u=3;
d=0.06;
a=c.*l.*(u.*d).^n./(d.*ni.^n);

function b=Biot(a,l,X)

b=a.*X./l;

function f=Fourier(l,ro,cp,t,x)
a=l./(ro.*cp);
f=a.*t./(x.^2);

function x_=Xbez(x1)
x_=x1./2;

function u=betan(Bi,n)

l=0;
x0=pi./3;
for i=1:n
    u(1,i)=fsolve(@(x) mi(Bi,x),x0+l.*pi);
    clc
    l=l+1;
end

function m=mi(Bi,x)

m=cot(x)-x./Bi;

function teta=rpc(Bi,x,Fo)

b=betan(Bi,20);
s=0;
for i=1:20
    s=s+An(Bi,i).*exp(-b(1,i).^2.*Fo);
end
teta=s;

function a=An(Bi,n)

b=betan(Bi,n);

a=2.*sin(b(1,n))./(b(1,n)+sin(b(1,n)).*cos(b(1,n)));

function f=rownanie(Bi1,Bi2,N,X,Fo,x)

Tp=396+273.15;
Tk=359+273.15;
T0=22+273.15;

f=tbezw(Tp,Tk,T0)-rpc(Bi1,2.*x/X,Fo).^2.*rpc(Bi2,2.*x./(N.*X),Fo);