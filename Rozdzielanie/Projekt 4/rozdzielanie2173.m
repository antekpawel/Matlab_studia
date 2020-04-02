%% Pawe³ Antkowiak
clear
close all
clc

%% dane
n1=4;
n2=4;
s=1e4/3600;
xs=1.8e-4;
m1=160;
m2=360;
d=3000/3600;
g=25/3600;


%% Obliczenia
w=s-d;
G=1e-4:1e-4:5*g;

for i=1:max(size(G))
    F=fsolve(@(x) dzban(m1,m2,xs,s,n1,n2,d,G(i),x),[0 0]);
    xd(i)=F(1);
    xw(i)=F(2);
    clc
end

plot1(G,xd,1,'b')
plot(G,xw,'r')

yd=m2*interp1(G,xd,g,'PCHIP')
yw=m1*interp1(G,xw,g,'PCHIP')




%% Funkcje

function F=dzban(m1,m2,xs,s,n1,n2,d,g,x)
w=s-d;
ag=((s/m1/g)^(n1+1)-s/m1/g)/((s/m1/g)^(n1+1)-1);
dg=((s/m2/g)^(n2+1)-s/m2/g)/((s/m2/g)^(n2+1)-1);
yw=(x(2)-x(1)*(1-dg))/(dg)*m2;
yd=(yw-ag*m1*xs)/(1-ag);
F(1)=g*(yd-yw)-s*(x(1)-xs);
F(2)=g*(yd-yw)-w*(x(1)-x(2));
end



