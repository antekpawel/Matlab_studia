function Projekt_2_zad1

%ZAD1

K=        273.15;
bar1=     1e5;
mmHg=     133.3224;
kPa=      1e3;
MPa=      1e6;
pa=       101325;
k=        1000;

%data

%zad1
zad1_p1=  15*kPa;
zad1_p2=  6*MPa;
zad1_p3=  200*kPa;
zad1_t=   450+K;

%import
[zad1_t_liq_sat,zad1_s_liq_sat,zad1_t_vap_sat,zad1_s_vap_sat]=WaterProperties;
[zad1_t_isobar_low,zad1_s_isobar_low,zad1_h_isobar_low]=WaterProperties_isop_015;
[zad1_t_isobar_high,zad1_s_isobar_high,zad1_h_isobar_high]=WaterProperties_isop_6;
[zad1_t_isobar_medium,zad1_s_isobar_medium,zad1_h_isobar_medium]=WaterProperties_isop_02;

zad1_enthropy_lowp_liq=0.75486;
zad1_enthropy_lowp_vap=8.0071;
zad1_enthropy_medp_vap=7.1269;
zad1_enthropy_medp_liq=1.5302;

zad1_enthalpy_3=2784.6;
zad1_enthalpy_4=2448;
zad1_enthalpy_1=225.93;
zad1_enthalpy_2=interp1(zad1_s_isobar_high,zad1_h_isobar_high,zad1_enthropy_lowp_liq,'spline');
zad1_enthalpy_2p=1213.9;

zad1_tem_saturation_lowp= 327.12;
zad1_tem_saturation_highp=548.73;
zad1_tem_saturation_medp=393.36;

zad1_enthalpy_lowp_vap=2598.3;
zad1_enthalpy_lowp_liq=225.94;

zad1_enthalpy_medp_vap=2706.2;
zad1_enthalpy_medp_liq=504.70;


%computing
zad1_temperature_4=interp1(zad1_s_isobar_high,zad1_t_isobar_high,zad1_enthropy_lowp_liq);
zad1_enthropy_highp_overheated=interp1(zad1_t_isobar_high,zad1_s_isobar_high,zad1_t,'spline');
zad1_enthropy_medp_overheated=interp1(zad1_t_isobar_medium,zad1_s_isobar_medium,zad1_t,'spline');
zad1_tem_lowp_overheated=interp1(zad1_s_isobar_low,zad1_t_isobar_low,zad1_enthropy_medp_overheated,'spline');

zad1_x_dry_a=(5.8901-zad1_enthropy_lowp_liq)/(zad1_enthropy_lowp_vap-zad1_enthropy_lowp_liq);
zad1_enthalpy_a=zad1_x_dry_a.*zad1_enthalpy_lowp_vap+(1-zad1_x_dry_a).*zad1_enthalpy_medp_liq;
zad1_turbine_work=zad1_enthalpy_a-zad1_enthalpy_3;
zad1_pump_work=zad1_enthalpy_2-zad1_enthalpy_lowp_liq;
zad1_heat=zad1_enthalpy_3-zad1_enthalpy_1;
zad1_heat_liq=zad1_enthalpy_lowp_liq-zad1_enthalpy_a;
zad1_efficiency=zad1_turbine_work/zad1_heat;

zad1_x_dry=(zad1_enthropy_highp_overheated-zad1_enthropy_lowp_liq)/(zad1_enthropy_lowp_vap-zad1_enthropy_lowp_liq);
zad1_enthalpy_b_hp_interp=interp1(zad1_s_isobar_high,zad1_h_isobar_high,zad1_enthropy_highp_overheated,'spline');
zad1_enthalpy_twophase=zad1_enthalpy_lowp_vap*zad1_x_dry+(1-zad1_x_dry)*zad1_enthalpy_lowp_liq;
zad1_turbine_work_b=zad1_enthalpy_b_hp_interp-zad1_enthalpy_twophase;
zad1_heat_b=zad1_enthalpy_b_hp_interp-zad1_enthalpy_2;
zad1_heat_liq_b=zad1_enthalpy_1-zad1_enthalpy_twophase;
zad1_efficiency_b=zad1_turbine_work_b/zad1_heat_b;

zad1_x_dry_c=(zad1_enthropy_highp_overheated-zad1_enthropy_medp_liq)/(zad1_enthropy_medp_vap-zad1_enthropy_medp_liq);
zad1_enthalpy_twophase_c=zad1_enthalpy_medp_vap*zad1_x_dry_c+(1-zad1_x_dry_c)*zad1_enthalpy_medp_liq;
zad1_enthalpy_medp_oh=interp1(zad1_s_isobar_medium,zad1_h_isobar_medium,zad1_enthropy_medp_overheated,'spline');
zad1_enthalpy_lowp_oh=interp1(zad1_s_isobar_low,zad1_h_isobar_low,zad1_enthropy_medp_overheated,'spline');
zad1_turbine_work_c=(zad1_enthalpy_b_hp_interp-zad1_enthalpy_twophase_c)+(zad1_enthalpy_medp_oh-zad1_enthalpy_lowp_oh);
zad1_heat_c=(zad1_enthalpy_b_hp_interp-zad1_enthalpy_2)+(zad1_enthalpy_medp_oh-zad1_enthalpy_medp_vap);
zad1_heat_liq_c=zad1_enthalpy_1-zad1_enthalpy_lowp_oh;
zad1_efficiency_c=zad1_turbine_work_c/zad1_heat_c;


%ploting

figure(11)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor;
hold all
xlabel('Enthropy [J/(g*K)]');
ylabel('Temperature [K]');
title('Rankine cycle');

plot(zad1_s_liq_sat,zad1_t_liq_sat,'g');
plot(zad1_s_isobar_low,zad1_t_isobar_low,'Color',[1 0.6 0]);
plot(zad1_s_isobar_high,zad1_t_isobar_high,'r');
plot([3.0278 5.8901 5.8901 zad1_enthropy_lowp_liq zad1_enthropy_lowp_liq],[548.73 548.73 327.12 327.12 zad1_temperature_4],'k--','LineWidth',1.5);%first proces

legend('Saturation Line','Isobaric 15kPa','Isobaric 20MPa','Cycle','Location','Best');
plot(zad1_s_vap_sat,zad1_t_vap_sat,'g')

plot(5.8901,548.73,'kv','MarkerFaceColor','k');%first point
plot(5.8901,327.12,'k<','MarkerFaceColor','k');%second point
plot(zad1_enthropy_lowp_liq,327.12,'k^','MarkerFaceColor','k');%third point
plot(zad1_enthropy_lowp_liq,zad1_temperature_4,'kx');%fourth point
plot(3.0278,548.73,'k>','MarkerFaceColor','k');%fifth point


j=1;h=1;
for i=1:max(size(zad1_s_isobar_high))
    if zad1_s_isobar_high(i)>zad1_enthropy_lowp_liq && zad1_s_isobar_high(i)<3.0278
        zad1_line_s(j)=zad1_s_isobar_high(i);
        zad1_line_t(j)=zad1_t_isobar_high(i);
        j=j+1;
    end
end
plot(zad1_line_s,zad1_line_t,'k--','LineWidth',1.5);

saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad1 a all','epsc')

axis([0.75 0.76 326.4 327.6]);

saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad1 point','epsc')


figure(12)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor;
hold all;
xlabel('Enthropy [J/(g*K)]');
ylabel('Temperature [K]');
title('Rankine cycle');

plot(zad1_s_liq_sat,zad1_t_liq_sat,'g')
plot(zad1_s_isobar_low,zad1_t_isobar_low,'Color',[1 0.6 0]);
plot(zad1_s_isobar_high,zad1_t_isobar_high,'r');
plot([zad1_enthropy_highp_overheated zad1_enthropy_highp_overheated],[zad1_tem_saturation_lowp zad1_t],'k--','LineWidth',1.5);

legend('Saturation Line','Isobaric 15kPa','Isobaric 20MPa','Cycle','Location','Best');
plot(zad1_s_vap_sat,zad1_t_vap_sat,'g')

plot(zad1_enthropy_lowp_liq,zad1_temperature_4,'kx');
plot(3.0278,548.73,'k>','MarkerFaceColor','k');
plot(zad1_enthropy_lowp_liq,327.12,'k^','MarkerFaceColor','k');

j=1;
for i=1:max(size(zad1_s_isobar_high))
    if zad1_s_isobar_high(i)>zad1_enthropy_lowp_liq && zad1_s_isobar_high(i)<3.0278
        zad1_line_s(j)=zad1_s_isobar_high(i);
        zad1_line_t(j)=zad1_t_isobar_high(i);
        j=j+1;
    end
end
plot(zad1_line_s,zad1_line_t,'k--','LineWidth',1.5);


plot(zad1_enthropy_highp_overheated,zad1_t,'kv','MarkerFaceColor','k');
plot(zad1_enthropy_highp_overheated,zad1_tem_saturation_lowp,'k<','MarkerFaceColor','k');
plot([zad1_enthropy_highp_overheated zad1_enthropy_lowp_liq],[zad1_tem_saturation_lowp zad1_tem_saturation_lowp],'k--','LineWidth',1.5);

j=1;
for i=1:max(size(zad1_s_isobar_high))
    if zad1_s_isobar_high(i)>zad1_enthropy_lowp_liq && zad1_s_isobar_high(i)<zad1_enthropy_highp_overheated
        zad1_b_left_s(j)=zad1_s_isobar_high(i);
        zad1_b_left_t(j)=zad1_t_isobar_high(i);
        j=j+1;
    end
end
plot(zad1_b_left_s,zad1_b_left_t,'k--','LineWidth',1.5);
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad1 b','epsc')

figure(13)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor
hold all
xlabel('Enthropy [J/(g*K)]');
ylabel('Temperature [K]');
title('Rankine cycle');

plot(zad1_s_liq_sat,zad1_t_liq_sat,'g');
plot(zad1_s_isobar_low,zad1_t_isobar_low,'Color',[1 0.6 0]);
plot(zad1_s_isobar_high,zad1_t_isobar_high,'r');
plot(zad1_s_isobar_medium,zad1_t_isobar_medium,'Color',[153 101 21]/255);
plot([zad1_enthropy_highp_overheated zad1_enthropy_highp_overheated zad1_enthropy_medp_vap],[zad1_t zad1_tem_saturation_medp zad1_tem_saturation_medp],'k--','LineWidth',1.5);

legend('Saturation Line','Isobaric 15kPa','Isobaric 20MPa','Isobaric 200kPa','Cycle','Location','Best');
plot(zad1_s_vap_sat,zad1_t_vap_sat,'g')

plot(5.8901,548.73,'kv','MarkerFaceColor','k');%first point
plot(zad1_enthropy_lowp_liq,327.12,'k^','MarkerFaceColor','k');%third point
plot(zad1_enthropy_lowp_liq,zad1_temperature_4,'kx');
plot(3.0278,548.73,'k>','MarkerFaceColor','k');

j=1;
for i=1:max(size(zad1_s_isobar_high))
    if zad1_s_isobar_high(i)>zad1_enthropy_lowp_liq && zad1_s_isobar_high(i)<3.0278
        zad1_line_s(j)=zad1_s_isobar_high(i);
        zad1_line_t(j)=zad1_t_isobar_high(i);
        j=j+1;
    end
end
plot(zad1_line_s,zad1_line_t,'k--','LineWidth',1.5);

plot(zad1_enthropy_highp_overheated,zad1_t,'kv','MarkerFaceColor','k');
plot([zad1_enthropy_lowp_vap zad1_enthropy_lowp_liq],[zad1_tem_saturation_lowp zad1_tem_saturation_lowp],'k--','LineWidth',1.5);
plot(zad1_enthropy_highp_overheated,zad1_tem_saturation_medp,'k>','MarkerFaceColor','k');
plot(zad1_enthropy_medp_overheated,zad1_t,'kv','MarkerFaceColor','k');
plot(zad1_enthropy_medp_overheated,zad1_tem_lowp_overheated,'k^','MarkerFaceColor','k');
%text(zad1_enthropy_medp_overheated,zad1_tem_lowp_overheated,'\leftarrow 1');
plot([zad1_enthropy_medp_overheated zad1_enthropy_medp_overheated],[zad1_tem_lowp_overheated zad1_t],'k--','LineWidth',1.5);

j=1;
for i=1:max(size(zad1_s_isobar_high))
    if zad1_s_isobar_high(i)>zad1_enthropy_lowp_liq && zad1_s_isobar_high(i)<zad1_enthropy_highp_overheated
        zad1_b_left_s(j)=zad1_s_isobar_high(i);
        zad1_b_left_t(j)=zad1_t_isobar_high(i);
        j=j+1;
    end
end
plot(zad1_b_left_s,zad1_b_left_t,'k--','LineWidth',1.5);

j=1;
for i=1:max(size(zad1_s_isobar_medium))
    if zad1_s_isobar_medium(i)>zad1_enthropy_medp_vap && zad1_s_isobar_medium(i)<zad1_enthropy_medp_overheated
        zad1_c_med_s(j)=zad1_s_isobar_medium(i);
        zad1_c_med_t(j)=zad1_t_isobar_medium(i);
        j=j+1;
    end
end
plot(zad1_c_med_s,zad1_c_med_t,'k--','LineWidth',1.5);

j=1;
for i=1:max(size(zad1_s_isobar_low))
    if zad1_s_isobar_low(i)>zad1_enthropy_lowp_vap && zad1_s_isobar_low(i)<zad1_enthropy_medp_overheated
        zad1_c_left_s(j)=zad1_s_isobar_low(i);
        zad1_c_left_t(j)=zad1_t_isobar_low(i);
        j=j+1;
    end
end
plot(zad1_c_left_s,zad1_c_left_t,'k--','LineWidth',1.5);
saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad1 c','epsc')



%display
diary C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\DaneZad1.txt
diary on
disp('a');
disp('Stopieñ suchoœci: '+string(zad1_x_dry_a));
disp('Entalpia na lini: '+string(zad1_enthalpy_a));
disp('Praca pompowania: '+string(zad1_pump_work));
disp('Praca wykonana przez turbine: '+string(zad1_turbine_work));
disp('Ciep³o pobrane w kotle i przegrzewaczu: '+string(zad1_heat));
disp('Ciep³o oddane w skraplaczu: '+string(zad1_heat_liq));
disp('Sprawnoœæ: '+string(zad1_efficiency));
disp('b');
disp('Stopieñ suchoœci: '+string(zad1_x_dry));
disp('Praca wykonana przez turbine: '+string(zad1_turbine_work_b));
disp('Ciep³o pobrane w kotle i przegrzewaczu: '+string(zad1_heat_b));
disp('Ciep³o oddane w skraplaczu: '+string(zad1_heat_liq_b));
disp('Sprawnoœæ: '+string(zad1_efficiency_b));
disp('c');
disp('Stopieñ suchoœci: '+string(zad1_x_dry_c));
disp('Praca wykonana przez turbine: '+string(zad1_turbine_work_c));
disp('Ciep³o pobrane w kotle i przegrzewaczu: '+string(zad1_heat_c));
disp('Ciep³o oddane w skraplaczu: '+string(zad1_heat_liq_c));
disp('Sprawnoœæ: '+string(zad1_efficiency_c));

diary off












