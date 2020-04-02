% Projekt rozdzielanie
clear
clc

tau1=0.38*60;%s
tau2=4.8*60;%s
v1=3.83e-3;%m^3
v2=15.17e-3;%m^3

tau1p=0.79*60;%s
tau2p=3.47*60;%s
v1p=6.78e-3;%m^3
v2p=14.88e-3;%m^3

fl=3.4;%m^2
dpl=3.1*98066.5;%MPa
dplp=4.40*98066.5;%MPa

f=38.40;%m
tau0=48*60;%s
mil=620e-3;%Pa*s
mi=540e-3;%Pa*s
mim=330e-3;%Pa*s
dp=4.2*98066.5;%MPa


KC=fsolve(@(x) dzban1(x,tau1,tau2,v1,v2),[0 0]);
KCp=fsolve(@(x) dzban1(x,tau1p,tau2p,v1p,v2p),[0 0]);

kl=KC(1);
c=KC(2);
klp=KCp(1);
cp=KCp(2);

s=1-(log(klp/kl)/log(dplp/dpl));


k=kl.*(dp./dpl).^(1-s).*(f/fl).^2.*mil./mi;
ck=c*(dpl/dp)^s*f/fl;

x0=[0.681,3000];

xaxis=1:1e4;
yaxis=-ck+sqrt((2.*ck).^2+4.*xaxis)./2;
figure(1)
clf
hold all
grid on
grid minor
plot(xaxis,yaxis,'r');
plot(-tau0,0,'ko')

ylabel('$$ V $$','Interpreter','latex')
xlabel('$$ \tau $$','Interpreter','latex')

x=fsolve(@(x) dzban2(x,k,ck,tau0),x0)
plot([-tau0 1.2*x(2)],[0 118],'k-');

x(2)+tau0
x(1)/(x(2)+tau0)
k/(2*(x(1)+ck))*mi/mim

function F=dzban1(x,tau1,tau2,v1,v2)
F(1)=v1^2+2*x(2)*v1-x(1)*tau1;
F(2)=v2^2+2*x(2)*v2-x(1)*tau2;
%x1=v
%x2=tau
end

function F=dzban2(x,k,c,tau0)
F(1)=x(1)/(x(2)+tau0)-k/(2*(x(1)+c));
F(2)=x(1)^2+2*c*x(1)-k*x(2);
%x1=v
%x2=tau
end

