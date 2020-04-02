function Projekt_2_zad2

%ZAD2

K=        273.15;
bar1=     1e5;
MPa=      1e6;

%data

%R134a
zad2_t=   55+K;
zad2_p1=  0.9*bar1/MPa;
zad2_p2=  22*bar1/MPa;
zad2_sm=  100/60;

%import
[t_sat_liq,p_sat_liq,s_sat_liq,h_sat_liq,t_sat_vap,p_sat_vap,s_sat_vap,h_sat_vap]=R134aSaturationLine;
[t_isop_009,s_isop_009,h_isop_009]=R134aisop_009;
[t_isop_22,s_isop_22,h_isop_22,cp_isop_22]=R134aisop_22;
[t_isop_044497,s_isop_044497,h_isop_044497]=R134aisop_044497;

s_sat_vap_lowp=     1.7499;
h_sat_vap_lowp=     381.18;
t_sat_lowp=         244.52;
h_sat_liq_lowp=     162.54;
s_sat_liq_lowp=     0.85576;
s_sat_vap_highp=    1.6941;
h_sat_vap_highp=    428.84;
t_sat_highp=        344.88;
s_sat_liq_highp=    1.3417;
h_sat_liq_highp=    307.30;

[t_isop_water,cp_isop_water]=WaterProperties_isop_01;
s_sat_vap_medp=    1.7212;
h_sat_vap_medp=    405.50;
t_sat_medp=        285.28;
s_sat_liq_medp=    1.0587;
h_sat_liq_medp=    216.51;

%computing
t_after_decompress=interp1(s_isop_22,t_isop_22,s_sat_vap_lowp,'spline');
h_after_compress=interp1(s_isop_22,h_isop_22,s_sat_vap_lowp,'spline');
q_liq=h_sat_vap_lowp-h_sat_liq_highp;
q_vap=h_after_compress-h_sat_liq_highp;
l_compress_ideal=h_after_compress-h_sat_vap_lowp;
teretical_coiffiency=q_vap/l_compress_ideal;
l_compress_real=l_compress_ideal/0.85;
real_coiffiency=q_vap/l_compress_real;

q_water_kg=integral(@(t) interp1(t_isop_water,cp_isop_water,t,'spline'),293.15,zad2_t);
q_water_all=q_water_kg.*zad2_sm;
q_water_all_with_loss=q_water_all/0.85;
%q_r134a_kg=integral(@(t) interp1(t_isop_22,cp_isop_22,t,'spline'),t_sat_highp,zad2_t);
q_r134a_kg=h_sat_liq_highp-h_after_compress;
mass_flow_r134a=abs(q_water_all_with_loss/q_r134a_kg);
compressor_power=l_compress_real*mass_flow_r134a;

p_medium_avaragegeo=sqrt(zad2_p2*zad2_p1);
h_mixer_med=interp1(s_isop_044497,h_isop_044497,s_sat_vap_lowp,'spline');
h_mixer_high=interp1(s_isop_22,h_isop_22,s_sat_vap_medp,'spline');
l_compress_ideal_b=h_after_compress-h_sat_vap_medp+(h_mixer_med-h_sat_vap_lowp);
q_r134a_kg_b=h_sat_liq_highp-h_mixer_high;
l_compress_real_b=l_compress_ideal_b/0.85;
teretical_coiffiency_b=-q_r134a_kg_b/l_compress_ideal_b;
real_coiffiency_b=-q_r134a_kg_b/l_compress_ideal_b.*0.85;
mass_flow_r134a_b=abs(q_water_all_with_loss/q_r134a_kg_b);
compressor_power_b=l_compress_real_b*mass_flow_r134a_b;

teretical_coiffiency_Mach=(h_mixer_high-h_sat_liq_highp)/(h_mixer_med-h_sat_vap_lowp+(h_mixer_med-h_sat_liq_medp)/(h_sat_vap_medp-h_sat_liq_highp)*(h_mixer_high-h_sat_vap_medp));

%plotting
figure(21)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
grid on;
grid minor
hold all
xlabel('Enthropy [kJ/(kg*K)]');
ylabel('Temperature [K]');
title('Heat pump');

plot(s_sat_liq,t_sat_liq,'b');
plot(s_isop_009,t_isop_009,'r');
plot(s_isop_22,t_isop_22,'Color',[1 0.6 0]);

plot(s_sat_vap,t_sat_vap,'b');
plot(s_sat_vap_lowp,t_after_decompress,'^k','MarkerFaceColor','k');
plot([s_sat_vap_lowp s_sat_vap_lowp],[t_sat_lowp t_after_decompress],'k--','LineWidth',1.5);
plot(s_sat_vap_lowp,t_sat_lowp,'^k','MarkerFaceColor','k');
plot(s_sat_liq_highp,t_sat_highp,'vk','MarkerFaceColor','k');
plot([s_sat_liq_highp s_sat_vap_highp],[t_sat_highp t_sat_highp],'k--','LineWidth',1.5);

saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad2 TS','epsc')

figure(22)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1])
semilogy(h_sat_liq,p_sat_liq,'b');
grid on;
hold all
axis([50 500 0 100]);
xlabel('Enthalpy [kJ/kg]');
ylabel('Pressure [MPa]');
title('Heat pump');

semilogy([h_sat_liq_highp h_after_compress],[zad2_p2 zad2_p2],'Color',[1 0.6 0]);
semilogy([h_sat_vap_lowp h_sat_liq_lowp],[zad2_p1 zad2_p1],'r');
semilogy([h_sat_liq_highp h_sat_liq_highp],[zad2_p2 zad2_p1],'Color',[153 50 204]/255);
semilogy([h_sat_liq_highp h_sat_liq_highp],[zad2_p2 zad2_p1],'k--','LineWidth',1.5);

legend('Saturation Line','Isobaric 2.2 MPa','Isobaric 0.09 MPa','Isoenthalpy '+string(h_sat_liq_highp)+' kJ/kg','Cycle','Location','best');
semilogy(h_sat_vap,p_sat_vap,'b');

semilogy(h_after_compress,zad2_p2,'k<','MarkerFaceColor','k');
semilogy([h_sat_liq_highp h_after_compress],[zad2_p2 zad2_p2],'k--','LineWidth',1.5);
semilogy(h_sat_liq_highp,zad2_p2,'kv','MarkerFaceColor','k');
semilogy([h_sat_vap_lowp h_sat_liq_highp],[zad2_p1 zad2_p1],'k--','LineWidth',1.5);
semilogy(h_sat_liq_highp,zad2_p1,'k>','MarkerFaceColor','k');
semilogy(h_sat_vap_lowp,zad2_p1,'kv','MarkerFaceColor','k');
%semilogy([h_sat_vap_lowp h_after_compress],[zad2_p1 zad2_p2],'k--','LineWidth',1.5);
semilogy([h_sat_vap_lowp h_mixer_med h_after_compress],[zad2_p1 p_medium_avaragegeo zad2_p2],'k--','LineWidth',1.5);
text(h_mixer_med,p_medium_avaragegeo,'\leftarrow S=const');

saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad2 PH','epsc')

figure(23)
clf;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
semilogy(h_sat_liq,p_sat_liq,'b');
grid on;
hold all
axis([50 500 0 100]);
xlabel('Enthalpy [kJ/kg]');
ylabel('Pressure [MPa]');
title('Heat pump');

semilogy([h_sat_liq_highp h_after_compress],[zad2_p2 zad2_p2],'Color',[1 0.6 0]);
semilogy([h_sat_vap_lowp h_sat_liq_lowp],[zad2_p1 zad2_p1],'r');
semilogy([h_sat_liq_medp h_mixer_med],[p_medium_avaragegeo p_medium_avaragegeo],'g');
semilogy([h_sat_liq_highp h_sat_liq_highp],[zad2_p2 p_medium_avaragegeo],'Color',[153 50 204]/255);
semilogy([h_sat_liq_medp h_sat_liq_medp],[zad2_p1 p_medium_avaragegeo],'Color',[218 165 32]/255);
semilogy([h_sat_vap_lowp h_sat_liq_highp],[zad2_p1 zad2_p1],'k--','LineWidth',1.5);

legend('Saturation Line','Isobaric 2.2 MPa','Isobaric 0.09 MPa','Isobaric '+string(p_medium_avaragegeo)+' MPa','Isoenthalpy '+string(h_sat_liq_highp)+' kJ/kg','Isoenthalpy '+string(h_sat_liq_medp)+' kJ/kg','Cycle','Location','best');
semilogy(h_sat_vap,p_sat_vap,'b');

saveas(gcf,'C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\zad2 2PH','epsc')


%display
diary C:\Users\Pawe³\Desktop\Pulpit\Nauka\IIrok\Termodynamika2\DaneZad2.txt
diary on
disp('a');
disp('Ciep³o oddane w skraplaczu: '+string(q_liq));
disp('Ciep³o oddane w parowniku: '+string(q_vap));
disp('Praca potrzebna na sprê¿enie (idealna): '+string(l_compress_ideal));
disp('Teoretyczny wspolczynnik wydajnosci: '+string(teretical_coiffiency));
disp('Praca potrzebna na sprê¿enie (rzeczywista): '+string(l_compress_real));
disp('Rzeczywisty wspolczynnik wydajnosci: '+string(real_coiffiency));
disp('Ciep³o potrzebne do ogrzania kg wody: '+string(q_water_kg));
disp('Ciep³o potrzebne do ogrzania wody: '+string(q_water_all));
disp('Ciep³o potrzebne do ogrzania wody ze stratami: '+string(q_water_all_with_loss));
disp('Ciep³o oddawane przez R134a: '+string(q_r134a_kg));
disp('Strumieñ masy czynnika roboczego: '+string(mass_flow_r134a));
disp('Moc sprê¿arki: '+string(compressor_power));

disp('b');
disp('Cisnienie w obszarze posrednim: '+string(p_medium_avaragegeo));
disp('Praca potrzebna na sprê¿enie (idealna): '+string(l_compress_ideal_b));
disp('Teoretyczny wspolczynnik wydajnosci: '+string(teretical_coiffiency_b));
disp('Praca potrzebna na sprê¿enie (rzeczywista): '+string(l_compress_real_b));
disp('Rzeczywisty wspolczynnik wydajnosci: '+string(real_coiffiency_b));
disp('Ciep³o oddawane przez R134a: '+string(q_r134a_kg_b));
disp('Teoretyczny wspolczynnik wydajnosci(wzór Machniewskiego): '+string(teretical_coiffiency_Mach));
disp('Strumieñ masy czynnika roboczego: '+string(mass_flow_r134a_b));
disp('Moc sprê¿arki: '+string(compressor_power_b));


diary off






















