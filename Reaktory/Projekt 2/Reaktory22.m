% Pawe³ Antkowiak
clear
close all
clc

beta=1.5;
t=[100 300 900];
T=[13.1 22 39.3]+273.15;
alfa=[0.2 0.52 0.87];
ca0=2.8;

x=0.47;
% x=1;
k2=alfa./(ca0.*(1-alfa).*(beta-alfa).*t);

y=log(k2);
xwykres=1./T;

figure(1)
hold all
grid on
grid minor
xlabel('$$ 1\over T $$','Interpreter','latex')
ylabel('$$ lnk_2 $$','Interpreter','latex')
plot(xwykres,y,'bo')

dzban=polyfit(xwykres,y,1);
plot([3e-3 3.6e-3],[3e-3.*dzban(1)+dzban(2) 3.6e-3.*dzban(1)+dzban(2)],'r');
tx=texlabel('a='+string(dzban(1)));
tx2=texlabel(' b='+string(dzban(2)));
text(3.2e-3,-5,tx)
text(3.2e-3,-5.25,tx2)

k20=exp(dzban(2));
ea=-dzban(1).*8.314;

eap=ea.*x;
k2p=k20.*exp(-eap./8.314./T);
dz=1+beta+1./k2p./ca0./t;
for i=1:3
f(i)=fsolve(@(x) x^2-dz(i).*x+beta,0)
end






