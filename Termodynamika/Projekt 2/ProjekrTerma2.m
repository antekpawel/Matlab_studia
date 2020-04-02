%Pawe³ Antkowiak
%Projekt 2
%Termodynamika

close all;
clear;
clc;
delete('C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\DaneZad1.txt')
delete('C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\DaneZad2.txt')
delete('C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\DaneZad3.txt')
delete('C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\DaneZad4.txt')

%Const
K=        273.15;
bar1=     1e5;
mmHg=     133.3224;
kPa=      1e3;
MPa=      1e6;
pa=       101325;
k=        1000;

%zad1
zad1_p1=  15*kPa;
zad1_p2=  6*MPa;
zad1_p3=  200*kPa;
zad1_t=   450+K;

%zad2
%R134a
zad2_t=   55+K;
zad2_p1=  0.9*bar1;
zad2_p2=  22*bar1;
zad2_sm=  100/60;
zad2_t_wl=20+K;

%zad3
%metan
zad3_t=   10+K;
zad3_p=   200*bar1;

%zad4
zad4_p_d= 4*bar1;
zad4_x1=  0.6;
zad4_p1=  1.25*bar1;
zad4_tm=  25+K;

zad4_tab_from_dr=[
0	    0       76.70   25	0        0
0.0298	0.0823	74.9	0	0.0934	-0.127
0.0615	0.1555	73.1	0	0.1425	-0.180
0.1106	0.266	70.3	0	0.2838	-0.276
0.1435	0.3325	68.6	0	0.4001	-0.307
0.2585	0.495	63.8	0	0.5393	-0.311
0.3908	0.634	59.3	0	0.6693	-0.264
0.5318	0.747	55.3	0	0.8153	-0.179
0.663	0.829	52.3	0	1        0
0.7574	0.878	50.4	0	0        0
0.8604	0.932	48.5	0	0        0
1       1       46.3	0	0        0
];

zad4_ccl4_=[6.984733	1277.603	234.3475	50	150	-752.70	8.9661      -3.04E-02	3.45E-05		250	388	1.2756E+02];
zad4_cs2 =[7.182683	1297.190	255.2298	-22	 47	 54.124	0.07614                  0         0        250 400 7.3159E+01];

zad4_ccl4_abc_yaws=[7.217900	1303.790	254.3940	-111.57	278.85];
zad4_cs2_abc_yaws =[7.011440	1278.540	232.8880	-22.82	283.20];

%ZAD1
disp('Zad 1');
Projekt_2_zad1
disp(' ');

%ZAD2;
disp('Zad 2');
Projekt_2_zad2
disp(' ');

%ZAD3;
disp('Zad 3');
Projekt_2_zad3
disp(' ');

%ZAD4;
disp('Zad 4');
%Projekt_2_zad4
disp(' ');




