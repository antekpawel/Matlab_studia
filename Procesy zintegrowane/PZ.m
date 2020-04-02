%% Pawe³ Antkowiak
%Procesy zintegrowane

clear
clc
close all
K=273.15;
g=9.81;
c
% xlswrite('C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Tabele\Tab1',proba);
% saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Wykresy\Wyk1','emf')

%% Dane projektowe

k2=2e-6;
ke=0.5;
vr=0.2;
qd=2;
c0=40;
ca=40;
cb=ca;


%% Obliczenia

t=0:1e3:3e4-1;
licznik=(exp(k2.*2.*c0.*t./sqrt(ke))-1).*c0;
mianownik=(ke.^(-1/2)-1)+exp(k2.*2.*c0.*t./sqrt(ke)).*(ke.^(-1/2)+1);
cd=licznik./mianownik;

figure(1)
clf
hold all
grid on
grid minor
xlabel('$$ Czas [s] $$','Interpreter','latex')
ylabel('$$ Stezenie \left[ {mol \over m^3} \right] $$','Interpreter','latex')
title('Stezenie w funkcji czasu','Interpreter','latex')


plot(t,cd,'-b')
plot(t,c0-cd,'-r')

legend('c_C=c_D','c_A=c_B','Location','Best');
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Wyk1','emf')

tab1(:,1)=t;
tab1(:,2)=cd;
tab1(:,3)=c0-cd;

xlswrite('C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Tab1',tab1,'B2');

figure(2)
clf
hold all
grid on
grid minor
xlabel('$$ Czas [s] $$','Interpreter','latex')
ylabel('$$ Stopien przereagowania [-] $$','Interpreter','latex')
title('Stopien przereagowania w funkcji czasu','Interpreter','latex')

alfa=cd./c0;

plot(t,alfa,'r')

saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Wyk2','emf')


licznik1=(exp(k2.*2.*c0.*9e5./sqrt(ke))-1).*c0;
mianownik1=(ke.^(-1/2)-1)+exp(k2.*2.*c0.*9e5./sqrt(ke)).*(ke.^(-1/2)+1);
cd1=licznik1./mianownik1;

alfak=cd1/c0;

wydajnosc=cd1/c0;
tab3(:,1)=t;
tab3(:,2)=alfa;
%% Z adsorbentem

t=0:6.7e3:2e5-1;
ca=(k2.*t+1./c0).^-1;

figure(3)
clf
hold all
grid on
grid minor
xlabel('$$ Czas [s] $$','Interpreter','latex')
ylabel('$$ Stezenie \left[ {mol \over m^3} \right] $$','Interpreter','latex')
title('Stezenie w funkcji czasu','Interpreter','latex')

plot(t,c0-ca,'b');
plot(t,ca,'r');
legend('c_C=c_D','c_A=c_B','Location','Best');

saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Wyk3','emf')

figure(4)
clf
hold all
grid on
grid minor
xlabel('$$ Czas [s] $$','Interpreter','latex')
ylabel('$$ Stopien przereagowania [-] $$','Interpreter','latex')
title('Stopien przereagowania w funkcji czasu','Interpreter','latex')

alfa=(c0-ca)./c0;

plot(t,alfa,'r')
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Wyk4','emf')

cak=(k2.*Inf+1./c0).^-1;

n=vr*c0;
m=n/qd;

tab2(:,1)=t;
tab2(:,2)=c0-ca;
tab2(:,3)=ca;

xlswrite('C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Tab2',tab2,'B2');

tab4(:,1)=t;
tab4(:,2)=alfa;