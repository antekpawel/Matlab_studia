%% Pawe³ Antkowiak
%Projekt 4 
%Aparatura

clear
clc
close all
%% Obliczenia

vc=9.833;
l=3.7;
lk=0.5;

mdk=30.5;
mkk=116.5;
dwk=0.6;
dzk=0.614;

dw=1.6;
dz=1.634;
mdd=286;
mkd=476;

mwd=pi*(dz^2-dw^2)/4*7850*l
mwk=pi*(dzk^2-dwk^2)/4*7850*lk
mk=mdk+2*mkk+mwk
mpus=4*mkd+2*mdd+mwd+mwk
mpel=mpus+(995*0.35+720*99.65)*vc








