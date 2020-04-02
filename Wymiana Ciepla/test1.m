function [u1,u2,time,s,x3]=test1
0.09;9;359;396;22;
0.09;9;359;396;22;

Bi1=1.3324;
Bi2=11.9916;

u1=betan(Bi1,19);
u2=betan(Bi2,19);

time=fsolve(@(t) rownanie3(u1,u2,t),35)
x=linspace(-1,1,1e3);
x1=linspace(0,0.42849,1e3);
figure(1);
plot(x1,rownanie5(u1,u2,x),'b-')
xlabel('Odleg³oœæ [m]')
ylabel('Temperatura [K]')
title('Rozk³ad temperatury wzd³u¿ przek¹tnej')
grid on
axis([0,0.5,0,700]);

x3=-1:0.1:1;
s=rownanie5(u1,u2,x3);
x3=linspace(0,0.42849,22);


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

function te=rownanie1(u,t)
te=0;
a=4.444e-7;
X=0.09;


for i=1:19
    te=te+2.*(sin(u(1,i))./(u(1,i)+sin(u(1,i)).*cos(u(1,i)))).*exp(-(u(1,i).^2).*(a.*t./((X./2).^2)));
end

function te=rownanie2(u,t)
te=0;
a=4.444e-7;
X=0.09;

for i=1:19
    te=te+2.*((sin(u(1,i)))./(u(1,i)+sin(u(1,i)).*cos(u(1,i)))).*exp(-(u(1,i).^2).*(a.*t./((7*X./2).^2)));
end

function x=rownanie3(u1,u2,t)
teta=0.36328;

x=teta-(rownanie1(u1,t).^2).*rownanie2(u2,t);

function x=rownanie4(u1,u2,t)

x=(rownanie1(u1,t).^2).*rownanie2(u2,t);

function xe=rownanie5(u1,u2,x)

xe=(rownanie6(u1,x).^2).*rownanie7(u2,x);

xe=xe.*(22-396)+396;
xe=xe+273.15;


function te=rownanie6(u,x)
te=0;
a=4.444e-7;
X=0.09;
t=777.2503;

for i=1:19
    te=te+2.*((sin(u(1,i)).*cos(u(1,i).*x))./(u(1,i)+sin(u(1,i)).*cos(u(1,i)))).*exp(-(u(1,i).^2).*(a.*t./((X./2).^2)));
end

function te=rownanie7(u,x)

te=0;
a=4.444e-7;
X=0.09;
t=777.2503;

for i=1:19
    te=te+2.*((sin(u(1,i)).*cos(u(1,i).*x))./(u(1,i)+sin(u(1,i)).*cos(u(1,i)))).*exp(-(u(1,i).^2).*(a.*t./((7.*X./2).^2)));
end