%% Pawe³ Antkowiak
%Projekt 2 zad1
%Kinetyka

clear
clc
close all
%% Dane
a=3e-2;%[m]
b=5e-2;%[m]
c=9e-2;%[m]
wp=0.32;%[-]
Dab=1e-9;%[m2/s]
wr=0.05;
wk=0.27;
Bi=50;




%% Obliczenia
Cp=wp/(1+wp);
Ck=wk/(1+wk);
Cr=wr/(1+wr);

mi=Mip(Bi,10);
an=Anp(Bi,mi);

t=10;
y11=sum(an.*exp(-mi.^2.*Dab.*t./(b/2).^2));
y12=sum(an.*exp(-mi.^2.*Dab.*t./(c/2).^2));
y1k=y11.*y12;

C=y1k.*(Cp-Cr)+Cr;
W=C./(1-C)
%fsolve(@(x) dzban(x),1e5)










