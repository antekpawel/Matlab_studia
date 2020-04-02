%% Pawe³ Antkowiak
%Reaktory mgr. Projekt 1

clear
close all
clc

K=273.15;

%% Dane projetowe

% ts=291+K;
% R=1.905e-2;
% l=9.144e-2;
% cp=0.8374e3;
% ro=1.6;
% ui=91440e-2/3600;
% rob=1443;
% dh=-69714;
% kr=0.9353/3600*100;
% Dr=4645e-4/3600;
% ca0=5e-5*1e-3*1e6;
% t0=316.2+K;

ts=291+K;
R=1.905;
l=9.144;
cp=0.8374;
ro=1.6e-3;
ui=91440;
rob=1443e-3;
dh=-69714;
kr=0.9353;
Dr=4645;
ca0=5e-5;
t0=316.2+K;

%% Obliczenia

dr=R/16;
dl=l/96;

kol=max(size(0:dl:l));
wi=max(size(0:dr:R));

T_roz=zeros(wi,kol)+t0;
f_roz=zeros(wi,kol);

T_roz(:,1)=t0;
T_roz(wi,:)=ts;


for i=2:wi-1
    for j=2:kol
        T_roz(1,j)=T_roz(1,j-1)+kr*dl/ro/cp/ui/(dr^2)*4*(T_roz(2,j-1)-T_roz(1,j-1))-szybkosc(T_roz(1,j-1),f_roz(1,j-1))*rob*ca0*dh*dl/ro/cp/ui;
        f_roz(1,j)=f_roz(1,j-1)+Dr*dl/ui/(dr^2)*4*(f_roz(2,j-1)-f_roz(1,j-1))-szybkosc(T_roz(1,j-1),f_roz(1,j-1))*rob*dl/ui;
        T_roz(i,j)=T_roz(i,j-1)+kr*dl/ro/cp/ui/(dr^2)*(1/(i-1)*(T_roz(i+1,j-1)-T_roz(i,j-1))+T_roz(i+1,j-1)-2*T_roz(i,j-1)-T_roz(i-1,j))-szybkosc(T_roz(i-1,j),f_roz(i-1,j))*rob*ca0*dh*dl/ro/cp/ui;
        f_roz(i,j)=f_roz(i,j-1)+Dr*dl/ui*(1/((i-1)*dr)*((f_roz(i+1,j-1)-f_roz(i,j-1))/dr)+((f_roz(i+1,j-1)-2*f_roz(i,j-1)+f_roz(i-1,j-1))/dr^2))+szybkosc(T_roz(i-1,j),f_roz(i-1,j))*rob*dl/ui;
    end
end

[x,y]=meshgrid(linspace(0,R,kol),linspace(0,l,wi));
figure(1)
surf(x,y,T_roz);
colormap(jet)
colorbar
%% Funkcje
function R=Rf(ra,ca0,rob)
R=-ra*ca0*rob;
end

function ra=szybkosc(T,f)
ra=1223040*exp(-5389/T)*(1-f);
end





















