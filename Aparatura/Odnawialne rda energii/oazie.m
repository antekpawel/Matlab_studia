%% Pawe³ Antkowiak
%Odnawialne i alternatywne zrodla energii

clear
clc
close all
K=273.15;
g=9.81;
d=24.*3600;
% xlswrite('C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Tabele\Tab1',proba);
% saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\Studia\III rok\Procesy zintegrowane\Wykresy\Wyk1','emf')

%% Dane projektowe

Qnom=60.*(0.8:0.1:1.2);%[m^3./d]
row=1200;%[kg./m^3]
ron=1100;%[kg./m^3]
xw=0.05;%[kgsm./kg]
xorg_w=0.68;%[kg./kgsm]
xn=0.06;%[kgsm./kg]
xorg_n=0.72;%[kg./kgsm]
w=0.5;%objetosciowego
beta=0.9;%[-]
HdoD=0.8;%[-]
R=0.6;%[-]
kw=9.5e-7;%[1./s]
kn=4.5e-7;%[1./s]
dhr=7e5;%[J./kg]
Tf=38;%[°C]
Tp=130;%[°C]
z=0.92;%[m^3./kg]
xch4=0.7;%objetosciowe
etael=0.4;%[-]
etat=0.97;%[-]
xi=0.1;%[-]
p=40;%[kW]
k1=4;%[W./m^2K]
k2=2;%[W./m^2K]
x=150;%[-]
ucyrk=0.8;%[m./s]
etapomp=0.6;%[-]
Mc=0.15;%[-]
N=2;%[-]
fiDSI=500;%[-]

%% Dane znalezione/dobrane

D_rz_cyrk=248e-3;


%% Przeliczenie na SI

Qnom=Qnom./d;
p=p.*1000;
Tf=Tf+K;
Tp=Tp+K;

%% Bilans masowy

A=w.*Qnom.*row.*xw.*xorg_w;
B=Qnom.*row;
C=kw.*row;
D1=(1-w).*Qnom.*ron.*xn.*xorg_n;
E=Qnom.*ron;
F=kn.*ron;
G=(1-R).*(D1+A);

Vr=1995.83;


bm5(1)=Vr;
bm5(2)=D1(5)./(F.*Vr+E(5));
bm5(3)=A(5)./(C.*Vr+B(5));
bm5(4)=R;


tab_bm(:,5)=bm5;
tab_bm(1,(1:4))=Vr;
tab_bm(4,5)=R;


xwyl_w=A./(C.*Vr+B);
xwyl_n=D1./(F.*Vr+E);
R=1-(B.*xwyl_w+E.*xwyl_n)./(D1+A);

tab_bm(2,:)=xwyl_n;
tab_bm(3,:)=xwyl_w;
tab_bm(4,:)=R;

WLOT_s=A+D1;
WYLOT_s=B.*xwyl_w+E.*xwyl_n;
%% Wymiary komory fermentacyjnej 

V_zbiornik=Vr./beta;
D=((4.*V_zbiornik)./(pi.*HdoD))^(1./3);
H=HdoD.*D;
h=D./3^(1./2)./2;

%% Strumien produkowanego biogazu

del_m_org=R.*WLOT_s;
V_biogazgazu=del_m_org.*z;

%% pompy 

Q_pomp=1.2.*Qnom.*x; %[m^3./h]
Q1pompy=Q_pomp./N;
d_max=((4.*Q1pompy)./(pi.*ucyrk)).^(1./2); %[m]
Fi=(10.*w.*Qnom.*xw+(1-w).*Qnom.*xn)./(w.*Qnom+(1-w).*Qnom);
mi_zw=Mi(Tf).*(1+2.5.*Fi+10.05.*Fi.^2+0.0027.*exp(16.6.*Fi));

Re=ucyrk.*D_rz_cyrk.*Ro(Tf)./mi_zw;
L=3.*H./4+2;

%% Bilans energetyczny

T_wlot=(10:5:25)+K;
T_otoczenia=(T_wlot-K).*2-15+K;

a1=pi.*D.*beta.*H;
l=D./sqrt(3);
a2=pi.*D.*(1-beta).*H+pi.*D./2.*l;

dp_rurociag=L.*Ro(Tf).*g;

for i=1:5
    for j=1:4
%         be_para(i,j)=m_para.*i_para.*etat;
        be_wylot(i,j)=WYLOT_s(i).*Cp(Tf).*Tf;
        be_wlot(i,j)=WLOT_s(i).*Cp(10+K).*T_wlot(j);
        be_reakcja(i,j)=dhr.*del_m_org(i);
%         be_mieszanie(i,j)=N_mieszania(i).*Mc;
        be_straty(i,j)=k1.*a1.*(Tf-T_otoczenia(j))+k2.*a2.*(Tf-T_otoczenia(j));
    end
end



ucyrkrzecz











