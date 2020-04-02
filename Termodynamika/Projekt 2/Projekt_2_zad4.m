function Projekt_2_zad4

%ZAD4

%Const
K=        273.15;
bar1=     1e5;
mmHg=     133.3224;
kPa=      1e3;
MPa=      1e6;
pa=       101325;
k=        1000;

%data

p_d= 4*bar1;
x1_1=  0.6;
p1=  1.25*bar1;
tm=  25+K;
t_ref=0+K;

%import

tab_from_dr=[
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

tab_from_dr=[
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
];

ccl4_=[6.984733	1277.603	234.3475	50	150	-752.70	8.9661      -3.04E-02	3.45E-05		250	388	1.2756E+02];
cs2 =[7.182683	1297.190	255.2298	-22	 47	 54.124	0.07614                  0         0        250 400 7.3159E+01];

ccl4_yaws=[7.217900	1303.790	254.3940	-111.57	278.85];
cs2_yaws =[7.011440	1278.540	232.8880	-22.82	283.20];

xchart=[0
0.0934
0.1425
0.2838
0.4001
0.5393
0.6693
0.8153
1];
HM=[1896.83
2027.86
2096.75
2295.06
2458.33
2653.79
2836.38
3041.47
3301.00];

cp_ccl4=[54.124	0.07614	0	0];
cp_cs2=[-752.70	8.9661	-3.04E-02	3.45E-05];


%calculation
x1=tab_from_dr(:,1);
x2=1-tab_from_dr(:,1);
y1=tab_from_dr(:,2);
y2=1-tab_from_dr(:,2);
t_exp=tab_from_dr(:,3);

p_sat_ccl4=anthonie(ccl4_yaws,t_exp);
p_sat_cs2=anthonie(cs2_yaws,t_exp);

gam1=pa.*y1./x1./p_sat_ccl4;
gam2=pa.*y2./x2./p_sat_cs2;

sq_gam1=sqrt(gam1);
sq_gam2=sqrt(gam2);

vanLaar_con=polyfit(sq_gam1,sq_gam2,1);
a=vanLaar_con(1);
b=vanLaar_con(2);

B=b.^2;
A=B./a.^2;

vanLaar_con1=polyfit(x1./x2,1./sq_gam1,1);
a1=vanLaar_con1(1);
b1=vanLaar_con1(2);

A1=exp(1./b1.^2);
B1=sqrt(A1)./b1;

vanLaar_con2=polyfit(x2./x1,1./sq_gam2,1);
a2=vanLaar_con2(1);
b2=vanLaar_con2(2);

A2=exp(1./b2.^2);
B2=sqrt(A2)./b2;

A_ava=(A1+A2)./2;
B_ava=(B1+B2)./2;
A_ava=0.17187;
B_ava=0.23874;

[a,b,~,res]=vanLaar_findAB(y1,x1,p_sat_ccl4,t_exp);

GRTvanLaar_mach=GibbsRT_vanLaar(x1,A_ava,B_ava);
GRTvanLaar_my=GibbsRT_vanLaar(x1,a,b);
GRT=GibbsRT(x1,gam1,gam2);

res_my=sum(abs(res));
res_mach=sum(abs(GRTvanLaar_mach-GRT));

gam1_vl_mach=vanLaar1_t(x1,A_ava,B_ava,t_exp);
gam1_vl_my=vanLaar1_t(x1,a,b,t_exp);

gam2_vl_mach=vanLaar2_t(x1,A_ava,B_ava,t_exp);
gam2_vl_my=vanLaar2_t(x1,a,b,t_exp);

y1_mach=gam1_vl_mach.*p_sat_ccl4.*x1./pa;
y1_my=gam1_vl_my.*p_sat_ccl4.*x1./pa;

t_mach=fsolve(@(t) pressure(gam1_vl_mach,x1,anthonie(ccl4_yaws,t),gam2_vl_mach,anthonie(cs2_yaws,t)),t_exp);
t_my=fsolve(@(t) pressure(gam1_vl_my,x1,anthonie(ccl4_yaws,t),gam2_vl_my,anthonie(cs2_yaws,t)),x1);

tw_mach=interp1(x1,t_mach,x1_1,'spline');
tw_my=interp1(x1,t_my,x1_1,'spline');

q_1_mach=integral(@(t) cp(cp_ccl4,t),tm,tw_mach+K);
q_2_mach=integral(@(t) cp(cp_cs2,t),tm,tw_mach+K);
q_mach=x1_1.*q_1_mach+(1-x1_1).*q_2_mach;

q_1_my=integral(@(t) cp(cp_ccl4,t),tm,tw_my+K);
q_2_my=integral(@(t) cp(cp_cs2,t),tm,tw_my+K);
q_my=x1_1.*q_1_my+(1-x1_1).*q_2_my;


T_D=dkrysian(x1_1,ccl4_yaws,cs2_yaws,x1,interp1(x1,gam1_vl_my,x1_1),interp1(x1,gam2_vl_my,x1_1))
 
 
 
%plotting
figure(41)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor;
hold all
xlabel('x_1');
ylabel('\gamma');
title('Gamma interpolation');
axis([0 1 1 1.5]);

plot(x1,gam1,'xr','MarkerSize',13,'LineWidth',2);
plot(x1,gam2,'xb','MarkerSize',13,'LineWidth',2);
plot((min(x1):0.001:max(x1)),interp1(x1,gam1,(min(x1):0.001:max(x1)),'spline'),'r');
plot((min(x1):0.001:max(x1)),interp1(x1,gam2,(min(x1):0.001:max(x1)),'spline'),'b');

saveas(gcf,'C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad4 loggam','epsc')

figure(42)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor;
hold all
xlabel('x_1');
ylabel('\DeltaG^E/RT');
title('Calculation value versus experimental');
axis([0 1 0 max(GRTvanLaar_mach).*1.5]);

plot(x1,GRT,'xr','MarkerSize',13,'LineWidth',2);
plot((min(x1):0.001:max(x1)),interp1(x1,GRTvanLaar_mach,(min(x1):0.001:max(x1)),'spline'),'r');
plot((min(x1):0.001:max(x1)),interp1(x1,GRTvanLaar_my,(min(x1):0.001:max(x1)),'spline'),'b');

legend('Experimental points','Function from MathCad','Function from MatLab','Location','Best');

saveas(gcf,'C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad4 GRT','epsc')


figure(43)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor;
hold all
xlabel('x_1');
ylabel('y_1');
title('Calculation value versus experimental');

plot(x1,y1,'xr','MarkerSize',13,'LineWidth',2);
plot((min(x1):0.001:max(x1)),interp1(x1,y1_mach,(min(x1):0.001:max(x1)),'spline'),'b','LineWidth',2);
plot((min(x1):0.001:max(x1)),interp1(x1,y1_my,(min(x1):0.001:max(x1)),'spline'),'g','LineWidth',1.5);

legend('Experimental points','Function from MathCad','Function from MatLab','Location','Best');
saveas(gcf,'C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad4 y1','epsc')

figure(44)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor;
hold all
xlabel('x_1');
ylabel('Temperature [K]');
title('Calculation value versus experimental');

plot(x1,t_exp,'xr','MarkerSize',13,'LineWidth',2);
plot((min(x1):0.001:max(x1)),interp1(x1,t_mach,(min(x1):0.001:max(x1)),'spline'),'b','LineWidth',2);
plot((min(x1):0.001:max(x1)),interp1(x1,t_my,(min(x1):0.001:max(x1)),'spline'),'g','LineWidth',1.5);
plot(y1,t_exp,'xr','MarkerSize',13,'LineWidth',2);
plot((min(y1):0.001:max(y1)),interp1(y1,t_mach,(min(y1):0.001:max(y1)),'spline'),'b','LineWidth',2);
plot((min(y1):0.001:max(y1)),interp1(y1,t_my,(min(y1):0.001:max(y1)),'spline'),'g','LineWidth',1.5);



legend('Experimental points','Function from MathCad','Function from MatLab','Location','Best');

saveas(gcf,'C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad4 tem','epsc')


figure(45)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor;
hold all
xlabel('x_1');
ylabel('Residuum');
title('Error chart');
axis([0 1 -0.05 0.05]);

plot(x1,res,'xg','MarkerSize',13,'LineWidth',2);
plot(x1,(GRTvanLaar_mach-GRT),'xb','MarkerSize',13,'LineWidth',2);
plot([0 1],[0 0],'k--');


saveas(gcf,'C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad4 res','epsc')

figure(46)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor;
hold all
xlabel('x_1');
ylabel('Enthalpy [kJ/kg]');
title('Basic isotermy');


plot(xchart,HM);
saveas(gcf,'C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad4 isoterm','epsc')



%display


diary C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\DaneZad4.txt
diary on
disp('Sta貫:');
disp('Sta豉 A: '+string(A_ava));
disp('Sta豉 B: '+string(B_ava));
disp('Suma residu闚 moich: '+string(res_my));
disp('Suma residu闚 Machu: '+string(res_mach));
disp('T wrzenia Machu: '+string(tw_mach));
disp('T wrzenia ja: '+string(tw_my));
disp('Cieplo Mach: '+string(q_mach));
disp('Cieplo ja: '+string(q_my));

diary off

close all
%function
function p=anthonie(c,t)
mmHg=     133.3224;
p=10.^(c(1)-c(2)./(t+c(3)))*mmHg;


function [a,b,error,res]=vanLaar_findAB(p0,t,p1,t1)

%Relationship beetween vapour pressure and temperature lnPv=A+B/(T+C)
y = (p0)./p1./t.*101325./exp(1./t1);

%Anthonie function
F = @(x,t) exp(x(2)./(1+t.*x(1)./(1-t)).^2);

%Predictions
x0 = [0.2 0.2];

%Main
[x,error,res] = lsqcurvefit(F,x0,t,y);

%Constant
a=x(1);
b=x(2);

function gam=vanLaar1_t(x,A,B,T)

gam=exp(B./T./(1+x.*A./(1-x)).^2);

function gam=vanLaar2_t(x,A,B,T)

gam=exp(B.*A./T./(A+(1-x)./x).^2);


function G=GibbsRT_vanLaar(x,A,B)
G=A.*B.*x.*(1-x)./((A-B).*x+B);

function G=GibbsRT(x,gam1,gam2)
G=x.*log(gam1)+(1-x).*log(gam2);


function p=pressure(gam1,x1,p1,gam2,p2)
p=gam1.*x1.*p1+gam2.*(1-x1).*p2-101325;

function cp=cp(tab,t)
cp=tab(1)+tab(2).*t+tab(3).*t.*t+tab(4).*t.*t.*t;

















function T_D=dkrysian(x1,tab1,tab2,x,Gamma01,Gamma02)
P=4e5;
Tmax=100;
Tmin=50;
PC_WODA=45.6*1e5;
TC_WODA=556.35;
PC_A=79*1e5 ;
TC_A=552;
PC=x1*PC_WODA+(1-x1)*PC_A;
TC=x1*TC_WODA+(1-x1)*TC_A;
PD=4e5;
eps=0.1;
T=80;
for i=1:1000
pr=PD/PC;
tr=(T+273.15)/TC;
PD1=anthonie(tab1,T);
PD2=anthonie(tab2,T);
pr1=PD1/PC_WODA;
tr1=(T+273.15)/TC_WODA;
pr2=PD2/PC_A;
tr2=(T+273.15)/TC_A;
f=exp(9*pr/128/tr*(1-6/tr.^2));
f1=exp(9*pr1/128/tr1*(1-6/tr1.^2));
f2=exp(9*pr2/128/tr2*(1-6/tr2.^2));
Y=(f1*x*Gamma01*PD1+f2*(1-x)*Gamma02*PD2)/(f*4*10^5);
if abs(Y-1)<eps
    break
end
if (Y-1)>0
    Tmax=T;
    T=(Tmax+Tmin)/2;
else
    Tmin=T;
    T=(Tmin+Tmax)/2;
end
end
T_D=T;










