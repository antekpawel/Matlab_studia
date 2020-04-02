%% Pawe³ Antkowiak
%Projekt Krystalizator
%Aparaty

clear
close all
clc
g=9.81;
K=273.15;
%% Dane projektowe
rozmiar_oczek=[
2.36E-03
1.65E-03
1.17E-03
8.33E-04
5.89E-04
4.17E-04
2.95E-04
];
x_mas=[
0
0.017
0.163
0.3
0.325
0.161
0.034
];
tau=100*60;
dkr_podane=2e-3;
kxs=0.7;
mi=7e-3;
roc=1440;
ros=1840;
delta=8e-5;
alfa=0.85;
M=32/1000;
Ms=53/1000;
x_xnas=0.009;
u=0.15;


%% Obliczenia wykresy

% czas przebywania
t=0:1e2:4e4;
F=1-exp(-t./tau);
figure(1)
hold all
plot(t,F)
title('$$F(t)=1-exp\left(-{t \over \tau}\right)$$','Interpreter','latex')
xlabel('Czas [s]');
ylabel('Gêstoœæ rozk³adu');
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Aparatura\Projekt 5 krystalizator\Wykresy\Fgestoscirozkladu','emf')

f=1./tau.*exp(-t./tau);
figure(2)
hold all
plot(t,f)
title('$$f(t)={1 \over \tau}exp\left(-{t \over \tau}\right)$$','Interpreter','latex')
xlabel('Czas [s]');
ylabel('Gêstoœæ rozk³adu');
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Aparatura\Projekt 5 krystalizator\Wykresy\frozkladu','emf')

% rozmiar i liczba krysztalow
dsr(1)=0;
for i=2:7
dsr(i)=(rozmiar_oczek(i-1)+rozmiar_oczek(i))/2;
end
dsr=dsr';
n=x_mas./ros./pi.*6./dsr.^3;
n(1)=0;
nsum=sum(n);
x_licz=n./nsum;

%generowanie tabeli
tab1(:,1)=rozmiar_oczek*1000;
tab1(:,2)=x_mas;
tab1(:,3)=dsr*1000;
tab1(:,4)=n;
tab1(:,5)=x_licz;
xlswrite('C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Aparatura\Projekt 5 krystalizator\Tabele\tab1',tab1);

%rozklad masowy
figure(3)
hold all
bar(dsr*1000,x_mas,'b')
title('Rozk³ad masowy w zale¿noœci od œrednicy')
xlabel('Œrednia œrednica kryszta³ów [mm]');
ylabel('U³amek masowy');
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Aparatura\Projekt 5 krystalizator\Wykresy\Rozk³ad masowy w zale¿noœci od œrednicy','emf')

%rozklad liczbowy
figure(4)
hold all
bar(dsr*1000,x_licz,'b')
title('Rozk³ad liczbowy w zale¿noœci od œrednicy')
xlabel('Œrednia œrednica kryszta³ów [mm]');
ylabel('U³amek liczbowy');
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Aparatura\Projekt 5 krystalizator\Wykresy\Rozk³ad liczbowy w zale¿noœci od œrednicy','emf')

% srednica srednia
dksr=sum(dsr.*x_licz);

%wspolczynnik masy
kxd=delta./dksr.*alfa.*(roc.*dksr.*u./mi).^0.6.*(mi./M./delta).^0.3;
Kxm=1/(1/kxs+1/kxd);

%szybkoœæ wzrostu
G=2.*Kxm.*M.*x_xnas./ros;

%rozklad krysztalów w produkcie

ddk=G*tau;
dprod=dsr+ddk;

dk=0:1e-5:8e-3;
fdk=1./G./tau.*exp(-dk./G./tau);

figure(5)
hold all
plot(dk,fdk)
title('$$f(\Delta d_k)={1 \over \tau \cdot G}exp\left(-{\Delta d_k \over \tau \cdot G}\right)$$','Interpreter','latex')
xlabel('$$\Delta d_k [m]$$','Interpreter','latex');
ylabel('Gêstoœæ rozk³adu');
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Aparatura\Projekt 5 krystalizator\Wykresy\frozkladuprodukt','emf')

%% Prawdopodobienstwa
%wiersze -> œrednice, kolumny -> sita
err=1;dksrpom=dksr;err2=1;
while err2>1e-5
while err>1e-5
dproba=dsr(7);
proba(7,7)=integral(@(x) 1./G./tau.*exp(-x./G./tau),0,rozmiar_oczek(6)-dproba);
proba(7,6)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(6)-dproba,rozmiar_oczek(5)-dproba);
proba(7,5)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(5)-dproba,rozmiar_oczek(4)-dproba);
proba(7,4)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(4)-dproba,rozmiar_oczek(3)-dproba);
proba(7,3)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(3)-dproba,rozmiar_oczek(2)-dproba);
proba(7,2)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(2)-dproba,rozmiar_oczek(1)-dproba);
proba(7,1)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(1)-dproba,Inf);

dproba=dsr(6);
proba(6,7)=0;
proba(6,6)=integral(@(x) 1./G./tau.*exp(-x./G./tau),0,rozmiar_oczek(5)-dproba);
proba(6,5)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(5)-dproba,rozmiar_oczek(4)-dproba);
proba(6,4)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(4)-dproba,rozmiar_oczek(3)-dproba);
proba(6,3)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(3)-dproba,rozmiar_oczek(2)-dproba);
proba(6,2)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(2)-dproba,rozmiar_oczek(1)-dproba);
proba(6,1)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(1)-dproba,Inf);

dproba=dsr(5);
proba(5,7)=0;
proba(5,6)=0;
proba(5,5)=integral(@(x) 1./G./tau.*exp(-x./G./tau),0,rozmiar_oczek(4)-dproba);
proba(5,4)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(4)-dproba,rozmiar_oczek(3)-dproba);
proba(5,3)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(3)-dproba,rozmiar_oczek(2)-dproba);
proba(5,2)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(2)-dproba,rozmiar_oczek(1)-dproba);
proba(5,1)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(1)-dproba,Inf);

dproba=dsr(4);
proba(4,7)=0;
proba(4,6)=0;
proba(4,5)=0;
proba(4,4)=integral(@(x) 1./G./tau.*exp(-x./G./tau),0,rozmiar_oczek(3)-dproba);
proba(4,3)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(3)-dproba,rozmiar_oczek(2)-dproba);
proba(4,2)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(2)-dproba,rozmiar_oczek(1)-dproba);
proba(4,1)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(1)-dproba,Inf);

dproba=dsr(3);
proba(3,7)=0;
proba(3,6)=0;
proba(3,5)=0;
proba(3,4)=0;
proba(3,3)=integral(@(x) 1./G./tau.*exp(-x./G./tau),0,rozmiar_oczek(2)-dproba);
proba(3,2)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(2)-dproba,rozmiar_oczek(1)-dproba);
proba(3,1)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(1)-dproba,Inf);

dproba=dsr(2);
proba(2,7)=0;
proba(2,6)=0;
proba(2,5)=0;
proba(2,4)=0;
proba(2,3)=0;
proba(2,2)=integral(@(x) 1./G./tau.*exp(-x./G./tau),0,rozmiar_oczek(1)-dproba);
proba(2,1)=integral(@(x) 1./G./tau.*exp(-x./G./tau),rozmiar_oczek(1)-dproba,Inf);

dproba=dsr(1);
proba(1,7)=0;
proba(1,6)=0;
proba(1,5)=0;
proba(1,4)=0;
proba(1,3)=0;
proba(1,2)=0;
proba(1,1)=integral(@(x) 1./G./tau.*exp(-x./G./tau),0,Inf);

for i=1:7
    for j=1:7
        frakcje1(j,i)=n(j).*proba(j,i);
    end
end
for i=1:7
    frakcje(i)=sum(frakcje1(:,i));
end
x_molowe_prod=frakcje./sum(frakcje);

dksr_prod=sum((x_molowe_prod'.*dsr));

err=abs(dksr_prod-dksrpom);
dksrpom=dksr_prod
kxd=delta./dksrpom.*alfa.*(roc.*dksrpom.*u./mi).^0.6.*(mi./M./delta).^0.3;
Kxm=1/(1/kxs+1/kxd);
G=2.*Kxm.*M.*x_xnas./ros;

end

pro_zaw=(dkr_podane-rozmiar_oczek(2))/(rozmiar_oczek(1)-rozmiar_oczek(2));
n_zaw=frakcje(6)*pro_zaw;
frakcje2(7)=frakcje(7);
frakcje2(6)=frakcje(6)+n_zaw;
frakcje2(5)=frakcje(5)+n(5);
frakcje2(4)=frakcje(4)+n(4);
frakcje2(3)=frakcje(3)+n(3);
frakcje2(2)=frakcje(2)+n(2);
frakcje2(1)=frakcje(1)+n(1);

x_molowe2_prod=frakcje2./sum(frakcje2);

dksr_prod=sum((x_molowe2_prod'.*dsr));
err2=abs(dksr_prod-dksrpom);
dksrpom=dksr_prod;

kxd=delta./dksrpom.*alfa.*(roc.*dksrpom.*u./mi).^0.6.*(mi./M./delta).^0.3;
Kxm=1/(1/kxs+1/kxd);
G=2.*Kxm.*M.*x_xnas./ros;
n=frakcje2;
end






xlswrite('C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Aparatura\Projekt 5 krystalizator\Tabele\tab2',proba);





figure(6)
hold all
bar(rozmiar_oczek*1000,x_molowe_prod,'r')
title('Rozk³ad liczbowy w zale¿noœci od œrednicy')
xlabel('Œrednia œrednica kryszta³ów [mm]');
ylabel('U³amek liczbowy');
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Aparatura\Projekt 5 krystalizator\Wykresy\rozklad liczbowy produktu re','emf')



% close all















