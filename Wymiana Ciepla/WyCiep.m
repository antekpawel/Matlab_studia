%Pawe³ Antkowiak
%Projekt 1
%Wymiana ciep³a

clear
clc
close all

%Data drom table
%zad1
%Cu
X=0.09;
%szk³o zwierciadlane
lam_mirror=[0.78,0.884];
dif_coeff_mirror=[0.392,0.4444].*1e-6;
lam_quartz=1.338;
dif_coeff_quartz=0.856.*1e-6;
lam_normal=0.744;
dif_coeff_normal=0.4444.*1e-6;

dif_coeff_choose=0.4444.*1e-6;
lam_choose=0.861;

N=9;
Tk=359;
Tp=396;
Up=3;
Lcalk=22;
T0=294.15;
C=0.16;
n=0.65;
%zad2
T1=346;
dw=0.14;
lam1=0.009;
ZB=1.7;
Z=1507;

dz=1.075*dw;
Tavarage=(T0+Tk)/2;
di=dz:0.001:0.5;

%Gut?
Tsurface=(T0+Tp)/2;
dsub=X;
precision_long=20;
precision_short=20;

%table from Gogó³
T_table=[50,60,70,80,90,100,120,140,160,180,200]+273.15;
lam_table=[2.83,2.90,2.97,3.05,3.13,3.21,3.34,3.49,3.64,3.78,3.93].*10^-2;
ni_table=[17.95,18.97,20.02,21.09,22.10,23.13,25.45,27.80,30.09,32.49,34.85].*10^-6;

%Computing

disp('Zadanie 1');

tem_less=(Tk-Tp)./(T0-Tp);
disp('Temperatura bezwymiarowa: '+string(tem_less));

pile_lenght=X.*N;

disp('D³ugoœæ stosu: '+string(pile_lenght)+' m');

lamp=interp1(T_table,lam_table,Tsurface);
nip=interp1(T_table,ni_table,Tsurface);

disp('Lepkoœæ dynamiczna powietrza przy sciance: '+string(nip)+' m^2/s'); 
disp('Wspó³czynnik przewodzenia ciep³a powietrza przy sciance: '+string(lamp)+' W/(m*K)');

alpha1=zad1_alpha(C,lamp,dsub,Up,nip,n);
disp('Wspó³czynnik wnikania ciep³a do p³ytek: '+string(alpha1)+' W/(m^2*K)');

Bi_short=Biot(alpha1,X,lam_choose);
Bi_long=Biot(alpha1,pile_lenght,lam_choose);

disp('Liczba Biota dla sciany krótkiej: '+string(Bi_short));
disp('Liczba Biota dla sciany d³ugiej:  '+string(Bi_long));

mi_np1=wyciep_mi(Bi_long,precision_long);
anp1=wyciep_anp(Bi_long,mi_np1);
disp('Wspó³czynnik mi.np i A.np dla liczby Biota sciany d³ugiej: ');
disp([mi_np1',anp1]);

mi_np2=wyciep_mi(Bi_short,precision_short);
anp2=wyciep_anp(Bi_short,mi_np2);
disp('Wspó³czynnik mi.np i A.np dla liczby Biota scian krótkich: ');
disp([mi_np2',anp2]);

time1=tem_less_iteration(mi_np1,mi_np2);
disp('Ten pierdolony czas: '+string(time1));

Fo1=dif_coeff_choose.*time1./(pile_lenght./2).^2;
Fo2=dif_coeff_choose.*time1./(X./2).^2;
disp('Liczba Fouriera dla sciany dlugiej : '+string(Fo1));
disp('Liczba Fouriera dla sciany krótkiej: '+string(Fo2));

vel=Lcalk./time1;
disp('Predkoœæ: '+string(vel));
tete=(Tavarage-Tp)/(T0-Tp);
disp('Druga teta: '+string(tete));

%zad2
disp(' ');
disp('Zadanie 2');

g=9.81665;

kizo=koszt_izo(di,Z,dz);
figure (1);
hold all
plot(di,kizo);

keks=koszt_e(di);
plot(di,keks);

kcal=keks+kizo;
plot(di,kcal,'g');

axis([dz 0.4 0 700]);

k_opt=min(kcal);

for i=1:max(size(kcal))
     if min(kcal)==kcal(i)
         di_opt=di(i);
     end
end

disp('Œrednica optymalna: '+string(di_opt)+' m');
disp('Jednostkowy koszt izolacji dla tej srednicy: '+string(k_opt));

plot ([0 di_opt di_opt],[k_opt k_opt 0],'k--');
scatter(di_opt,k_opt,'kx');

xlabel('Œrednica izolacji [m]');
ylabel('Jednostkowy koszt izolacji')
legend('Koszty inwestycyjne','Koszty eksploatacyjne','Koszt ca³kowity','Minimum kosztów')
title('Koszty w zale¿noœci od gruboœci izolacji');




rezkladtemperatury(mi_np1,anp1,X,mi_np2,anp2,pile_lenght,dif_coeff_choose,time1,T0,Tp);

























