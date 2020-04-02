function Projekt_2_zad3

%ZAD3

K=        273.15;
bar1=     1e5;
mmHg=     133.3224;
kPa=      1e3;
MPa=      1e6;
pa=       101325;
k=        1000;

%data 
%metan
zad3_t=   10+K;
zad3_p=   200*bar1;

%data import
[t_interp_isot,h_interp_isot,s_interp_isot]=MethaneInterpHIsotermic;
[t_isobaric1,h_isobaric1,s_isobaric1]=MethaneIsobaric_1;
[t_saturationl,h_staturationl,s_staturationl,t_saturationv,h_staturationv,s_staturationv]=MethaneSaturationLine;
[t_isobaric200,h_isobaric200,s_isobaric200]=MethaneIsobaric_200;
[s_it,h_it]=MethaneIsotermic_10

h_liquid_lowp=-0.55734;
h_vapor_lowp=510.56;
t_saturation_lowp=111.51;

s_liq_lowp=-0.0049665;
s_vap_lowp= 4.5787;

%computing

%Enthalpy
h_isoterm_highpress=interp1(t_isobaric200,h_isobaric200,zad3_t,'spline');
h_isoterm_lowpress=interp1(t_interp_isot,h_interp_isot,zad3_t,'spline');

w=(h_isoterm_lowpress-h_isoterm_highpress)/(h_isoterm_lowpress-h_liquid_lowp);

h_isoenthalpy=(w)*h_liquid_lowp+(1-w)*h_vapor_lowp;
t_interp_isoe=interp1(h_isobaric200,t_isobaric200,h_isoenthalpy);

j=1;
for i=1:max(size(t_isobaric200))
    if t_isobaric200(i)>t_interp_isoe && t_isobaric200(i)<zad3_t
        t_plot_line1(j)=t_isobaric200(i);
        h_plot_line1(j)=h_isobaric200(i);
        j=j+1;
    end
end
j=1;
for i=1:max(size(t_isobaric1))
    if t_isobaric1(i)>t_saturation_lowp && t_isobaric1(i)<zad3_t
        t_plot_line2(j)=t_isobaric1(i);
        h_plot_line2(j)=h_isobaric1(i);
        j=j+1;
    end
end

%enthropy

s_isoterm_highpress=interp1(t_isobaric200,s_isobaric200,zad3_t,'spline');
s_isoterm_lowpress=interp1(t_interp_isot,s_interp_isot,zad3_t,'spline');

%question

s_isoe_lowp=(w)*s_liq_lowp+(1-w)*s_vap_lowp;
s_isoe_highp=interp1(h_isobaric200,s_isobaric200,h_isoenthalpy);

l_lost_valve=zad3_t*(s_isoe_highp-s_isoe_lowp);
l_min_liq=(h_liquid_lowp-h_isoterm_lowpress)-zad3_t*(s_liq_lowp-s_isoterm_lowpress);
l_technical=zad3_t*(s_isoterm_lowpress-s_isoterm_highpress)-(h_isoterm_lowpress-h_isoterm_highpress);
l_all_work_lost=l_technical-w*l_min_liq;
q_heatexchanger_lost=l_all_work_lost+l_lost_valve;


%ploting

figure(31);
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor;
hold all;
axis([-250 1000 0 400]);
title('Linde cycle');
xlabel('Enthalpy [kJ/kg]');
ylabel('Temperature [K]');

plot(h_isobaric1,t_isobaric1,'g');
plot(h_staturationl,t_saturationl,'b');
plot(h_isobaric200,t_isobaric200,'Color',[255 102 102]/255);

plot([h_isoterm_highpress h_isoterm_lowpress],[zad3_t zad3_t],'Color',[15 125 18]/255);
plot([h_isoenthalpy h_isoenthalpy],[t_saturation_lowp t_interp_isoe],'Color',[1 0.6 0]);
plot(h_plot_line1,t_plot_line1,'k--','LineWidth',1.5);


legend('Isobaric 0.1MPa','Saturation line','Isobaric high 20MPa','Isotermic 283.15K','Isoenthalpy 392.5091 kJ/kg','Cycle','Location','Best');
plot(h_staturationv,t_saturationv,'b');

plot([h_isoenthalpy h_isoenthalpy],[t_saturation_lowp t_interp_isoe],'k--','LineWidth',1.5)
plot([h_isoterm_highpress h_isoterm_lowpress],[zad3_t zad3_t],'k--','LineWidth',1.5);
plot([h_vapor_lowp h_isoenthalpy],[t_saturation_lowp t_saturation_lowp],'k--','LineWidth',1.5);
plot(h_plot_line2,t_plot_line2,'k--','LineWidth',1.5);

saveas(gcf,'C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad3 TH','epsc')

figure(32);
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
hold all;
grid on;
grid minor;
axis([-2 10 0 400]);
title('Linde cycle');
xlabel('Enthropy [kJ/(kg*K)]');
ylabel('Temperature [K]');

plot(s_isobaric1,t_isobaric1,'g');
plot(s_staturationl,t_saturationl,'b');
plot(s_isobaric200,t_isobaric200,'Color',[255 102 102]/255);

plot([s_isoterm_highpress s_isoterm_lowpress],[zad3_t zad3_t],'Color',[15 125 18]/255);

h=area([s_isoe_lowp s_isoe_lowp s_isoe_highp s_isoe_highp],[0 zad3_t zad3_t 0]);
h(1).FaceColor =[0,191,255]/255;

plot([s_liq_lowp s_liq_lowp s_isoterm_lowpress],[t_saturation_lowp zad3_t zad3_t]);


legend('Isobaric 0.1MPa','Saturation line','Isobaric high 20MPa','Isotermic 283.15K','Valve loss of work','Location','SouthEast');

plot(s_staturationv,t_saturationv,'b');

plot(s_isobaric1,t_isobaric1,'g');
plot(s_staturationl,t_saturationl,'b');
plot(s_isobaric200,t_isobaric200,'Color',[255 102 102]/255);

saveas(gcf,'C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad3 TS','epsc')

figure(33);
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
hold all;
grid on;
grid minor;
title('Linde cycle');
xlabel('Enthropy [kJ/(kg*K)]');
ylabel('Enthalpy [kJ/kg]');
axis([-5 10 0 1000])

plot(s_isobaric1,h_isobaric1,'g');
plot(s_staturationl,h_staturationl,'b');
plot(s_isobaric200,h_isobaric200,'Color',[255 102 102]/255);
plot(s_it,h_it,'r')
legend('Isobaric 0.1MPa','Saturation line','Isobaric high 20MPa','Isotermic 283.15K','Location','SouthEast');



plot(s_staturationv,h_staturationv,'b');

saveas(gcf,'C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad3 SH','epsc')


%disp
diary C:\Users\Pawe許Desktop\Pulpit\Nauka\IIrok\Termodynamika2\DaneZad3.txt
diary on
disp('U豉mek skroplonego gazu: '+string(w));
disp('Temperatura przed zaworem: '+string(t_interp_isoe));
disp('Praca stracona na zaworze: '+string(l_lost_valve));
disp('Minimalna praca skraplania kilograma gazu: '+string(l_min_liq));
disp('Sumaryczne straty pracy: '+string(l_all_work_lost));
disp('Praca techniczna: '+string(l_technical));
disp('Staraty ciep豉 w wymienniku: '+string(q_heatexchanger_lost));
diary off


























