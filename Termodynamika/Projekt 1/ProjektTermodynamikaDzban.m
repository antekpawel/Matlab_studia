%Thermodynamics Project 1
%Pawe³ Antkowiak
%Wariant 31

close all
clear
clc

%Notes

%https://app.knovel.com/web/toc.v/cid:kpYCPDCECD/viewerType:toc/root_slug:yaws-critical-property/url_slug:table-47-specific-heat?b-q=material_or_substance_name%3A%22nitrogen%22%20AND%20specific_heat_mf%3A%5B%20*%20TO%20*%20%5D&b-facet-selected=item_type_nospace%3Aigraph&b-group-by=true&b-dsQuery=%3Fmn%3Dnitrogen%3Aand%3Apn%3Dspecific%20heat&b-search-type=tech-reference&b-sort-on=default
%zad4 problemy z interpolacja plaszczyzny
%dyfuzja metoda wilkego chenga

%Date

%Addytional constant

K=273.15;%Celcius to Kelwin
R=8.3144;% Gas constant
Hg=133.3224;%mmHg to Pascals
bar1=100000;%Bars to Pascals

%zad1
%gas A cyclohexane
%gas B ethane
%zmieniam temperature bo cykloheksa nie wrze
tem1_mix1=0+K;
tem1_mix2=120+K;
p1_mix2=115*bar1;
q1_str=7.3*1000;
m1=50./3600;
xa=0.75;
xb=0.25;
p1=1*bar1;
tem1=0+K;
dw1=20/1000;

%ethane CAS 74-84-0 M=30.07 Tmin=150K Cpmin=38.66 Tmax=1500K Cpmax=145.9
%zc=0.279
M1_eth=29;
tab1_moduls_eth=[
28.7167713794797
7.34582726095694E-03
-4.54759070492041E-05
1.16406082454094E-07
-1.22458456787057E-10
5.90449193986651E-14
-1.08748042093876E-17];
tc_eth=132.45;
pc_eth=37.74.*bar1;
om1_eth=0.0377;

%cyclohexane CAS 110-82-7 M=84.161 Tmin=150K Cpmin=68.471 Tmax=1500K
%Cpmax=368.815 zc=0.273 Temperatyra wrzenia 60°C
M1_cychex=44;
tab1_moduls_cychex=[
23.5061038797879
3.80655731822383E-02
7.40233050458845E-05
-2.22713256071084E-07
2.34374551216989E-10
-1.14647749938175E-13
2.16814726839333E-17];
tc_cychex=304.21;
pc_cychex=73.83.*bar1;
om1_cychex=0.2276;

%cyclohexane
tab1_moduls_vis_cychex=[6.77e-8;0.8367;36.7;0];%Perry

%Ethane
tab1_moduls_vis_eth=[2.5906e-7;0.67988;98.902;0];%Perry

data_from_carr_chart=5.1;

%zad2
%rotwór C ethanol
tem2_p=0+K;
x2_mol=0.5;
tem2_k=20+K;
m2=10;
M2_water=18.015/1000;
M2_ethanole=46.069/1000;

%ict_heat capacity 44951_03B str25
tem2_itc=[3,23,41]+K;
xmol=[0,2.02,4.16,6.46,8.91,11.5,14.4,20.7,28.1,37,47.7,61,77.9,100];
matrix_from_itc=[4.211,4.316,4.379,4.395,4.362,4.282,4.186,3.925,3.622,3.367,3.132,2.802,2.568,2.263
                 4.182,4.236,4.274,4.312,4.324,4.299,4.241,4.077,3.847,3.588,3.329,3.040,2.760,2.417
                 4.174,4.241,4.291,4.320,4.324,4.299,4.257,4.107,3.898,3.655,3.404,3.132,2.861,2.601]'.*1000;

%zad3
%N-metyloanilina
tem3_zak=60+K:1:200+K;
tem3=100+K;
tw3=195.5+K;
tc3=701.55;
pc3=52.*bar1;
acentric3=0.475;
M3=107.155/1000;

%CAS=100-61-8

%Nelson, O.A.; Wales, H., Vapor Pressures and Boiling Points of Mono- and Dimthylanilines and Mono- and Diethylanilines, J. Am. Chem. Soc., 1925, 47, 867-872.
%4.99409
%2226.576
%-22.456
%tmin323
%tmax472.76
%plot(tem3_zak,log(10.^(4.99409-(2226.576./(tem3_zak-22.456)))),'r');

%Yaws' Handbook of Thermodynamic and Physical Properties of Chemical Compounds
%P=10^(a+(b./t)+c.*log(t)+d.*t+e.*t.*t)
%-7.2448
%-2964.8
%8.9163
%-2.0616E-02
%1.0451E-05
%tmin216.15
%tmax701.55
yaws3=log(10.^(-7.2448+(-2964.8./tem3_zak)+8.9163.*log10(tem3_zak)+-2.0616E-02.*tem3_zak+1.0451E-05.*tem3_zak.*tem3_zak));


%Perry str 99 experimental points
p3_vapour_pressure=[1,5,10,20,40,60,100,200,400,760];
t3_vapour_temperature=[36,62.8,76.2,90.5,106,115.8,129.8,149.3,172,195.5]+K;

%Clausiusa-Clapeyrona from dr Machniewski
tc3_for_cc=701.55;
pc3_for_cc=52*bar1;

%zad4
%roztwór E (NH4)2SO4
c4_salt=2;%mol/kg
tem4=15+K;
M4=132.134/1000;

%zad5
%octan n-butylu
tem5_p=10+K;
tem5_k=50+K;

M5=116.16./1000;
associate5=2.6;

%Share element
element_c=                  14.8/1000000;
element_h=                  3.7/1000000;
element_o_duble=            7.4/1000000;
element_o_ester_methyl=     9.1/1000000;
element_o_ester_ethyl=      9.9/1000000;
element_o_ester_high=       11/1000000;

%Amount
amount_c=6;
amount_h=12;
amount_o_double=1;
amount_o_ester_high=1;

%Water viscosity 
%ict_viscosity 44951_01A str10
viscosity_tem5_p=13.097/10000;
viscosity_tem5_k=5.492/10000;

%zad 1

disp('Zadanie 1');
figure(1);
hold all;
title('Ciep³o w³aœciwe przy sta³ym ciœnieniu w zale¿noœci temperatury');
xlabel('Temperatura [K]');
ylabel('Ciep³o w³aœciwe przy sta³ym ciœnieniu [J/K]');

%Dane i obliczenie ciep³a w³aœciwego
%Yaws' Critical Property Data for Chemical Engineers and Chemists

tem1_range=tem1_mix1:0.1:tem1_mix2;

cp_eth=zad1_heat_gas(tem1_range,tab1_moduls_eth);
cp_cychex=zad1_heat_gas(tem1_range,tab1_moduls_cychex);

%Enthalpy
xma=xa.*M1_cychex/(xa.*M1_cychex+xb.*M1_eth);
xmb=xb.*M1_eth/(xa.*M1_cychex+xb.*M1_eth);

%Plotting
bar(tem1_range,cp_eth,'g');
bar(tem1_range,cp_cychex);

enthalpy_correction_pg=zad1_p_correction_pg(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb,om1_cychex,om1_eth);
enthalpy_sum=integral(@(t) zad1_heat_gas(t,tab1_moduls_cychex),tem1_mix1,tem1_mix2).*xma+integral(@(t) zad1_heat_gas(t,tab1_moduls_eth),tem1_mix1,tem1_mix2).*xmb+enthalpy_correction_pg;
enthalpy_sum_nocorect=integral(@(t) zad1_heat_gas(t,tab1_moduls_cychex),tem1_mix1,tem1_mix2).*xma+integral(@(t) zad1_heat_gas(t,tab1_moduls_eth),tem1_mix1,tem1_mix2).*xmb;

disp('Entalpia mieszaniny bez poprawki: '+string(enthalpy_sum_nocorect./1000)+' kJ/kg');
disp('Poprawka na entalpie: '+string(enthalpy_correction_pg./1000)+' kJ/kg');
disp('Ca³kowita entalpia: '+string(enthalpy_sum./1000)+' kJ/kg');

%Efficiency and power calculations

power1_compress_total=enthalpy_sum.*m1;
disp('Moc potrzebna do sprê¿enia: '+string(power1_compress_total));

spr=power1_compress_total/(power1_compress_total+q1_str);
disp('Sprawnoœæ '+string(spr));
Compressor_power=power1_compress_total+q1_str;
disp('Moc potrzebna do sprê¿enia '+string(Compressor_power)+'W');

%Kay's rule

%Check
check_kay_rule=[tc_cychex/tc_eth,pc_cychex/pc_eth];
disp('Sprawdzenie warunków regu³y Kaya: ');
disp('Cisnienia:   '+string(check_kay_rule(2)));
disp('Temperatury: '+string(check_kay_rule(1)));
if check_kay_rule(1)<=2 && check_kay_rule(1)>=0.5 && check_kay_rule(2)<=2 && check_kay_rule(2)>=0.5
    disp('Warunki regu³y Kaya s¹ spe³nione')
else
    disp('Warunki regu³y Kaya NIE s¹ spe³nione');
end

%Computing
tc_kay=tc_cychex.*xa+tc_eth.*xb;
pc_kay=pc_cychex.*xa+pc_eth.*xb;
disp('Warunki psuedokrytyczne wedlug regu³y Kaya: ');
disp('Cisnienie:   '+string(pc_kay)+'  Nasze warunki: '+string(p1_mix2));
disp('Temperatura: '+string(tc_kay)+'   Nasze warunki: '+string(tem1_mix2));

%Reducted
p_red=p1_mix2./pc_kay;
t_red=tem1_mix2./tc_kay;
disp('Parametry zredukowane:')
disp('Cisnienie:   '+string(p_red)+'   Temperatura: '+string(t_red));

%Communicat
if tc_kay<=tem1_mix2 && pc_kay<=p1_mix2
    disp('Warunki zgodne mieszanina znajduje siê w stanie pseudokrytycznym.');
else
    disp('Mieszanina nie znajdujê siê w stanie pseudokrytycznym');
end

%Dr Machniewski's function which fix roots from equation z=f(p,t,xa)
disp('Wspó³czynnik œciœliwoœci wyznaczony metod¹ Kaya z r. ');
%Redlich-Kwong
z_rk=zad1_redlich_kwong_2mix(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb);
disp('Redlicha-Kwonga: '+string(z_rk(3)));
%Penga-Robinsona
z_pr=zad1_peng_robinson_2mix(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb,om1_cychex,om1_eth);
disp('Penga-Robinsona: '+string(z_pr(3)));
%Van der Wals... please.
%Soave RKS
z_rks=zad1_soave_2mix(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb,om1_cychex,om1_eth);
disp('RKS:             '+string(z_rks(1)));

disp('Wspó³czynnik œciœliwoœci wyznaczony regu³¹ mieszania z r. ');
%Redlich-Kwong
z_rk_m=zad1_redlich_kwong_2mix_mixrule(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb);
disp('Redlicha-Kwonga: '+string(z_rk_m(3)));
%Penga-Robinsona
z_pr_m=zad1_peng_robinson_2mix_mixrule(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb,om1_cychex,om1_eth);
disp('Penga-Robinsona: '+string(z_pr_m(1)));
%Van der Wals... please.
%Soave RKS
z_rks_m=zad1_soave_2mix_mixrule(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb,om1_cychex,om1_eth);
disp('RKS:             '+string(z_rks_m(3)));

disp('Ze wzglêdu na warunki zadania wybieram metode Penga-Robinsona z metod¹ Kaya');

%Density
density1=(xa.*M1_cychex+xb.*M1_eth).*p1_mix2./1000./R./(tem1_mix2)./z_pr_m(3);
disp('Gêstoœæ mieszaniny wynosi: '+string(density1)+' kg/m3');

%Velocity
velocity1=m1./density1.*4./dw1./dw1./pi;
disp('Prêdkoœæ mieszaniny w przewodzie dw=20mm wynosi: '+string(velocity1)+'m/s');

%Viscosity 

Viscosity1_cychex=zad1_viscosity_gas_perry(tem1_mix2,tab1_moduls_vis_cychex);
disp('Lepkoœæ czystego cykloheksanu: '+string(Viscosity1_cychex)+'Pa/s');
Viscosity1_eth=zad1_viscosity_gas_perry(tem1_mix2,tab1_moduls_vis_eth);
disp('Lepkoœæ czystego etanu:        '+string(Viscosity1_eth)+'Pa/s');

%Viscosity for mixture
disp('Lepkoœæ mieszaniny wyznaczona metod¹:');
disp('Data odczytana z wykresu Carra: '+string(data_from_carr_chart));
%Herning-Ziperrer
Viscosity1_mix_hz=data_from_carr_chart.*zad1_viscosity_herning_zipperer(xa,xb,M1_cychex,M1_eth,tc_cychex,tc_eth,Viscosity1_cychex,Viscosity1_eth);
disp('Herninga-Zipperera: '+string(Viscosity1_mix_hz)+'Pa*s');
%Wilke
Viscosity1_mix_w=data_from_carr_chart.*zad1_viscosity_wilke(Viscosity1_cychex,Viscosity1_eth,xa,xb,M1_cychex,M1_eth);
disp('Wilkego:            '+string(Viscosity1_mix_w)+'Pa*s');

legend('Ethane','Cyclohexane');
%Enthropy
%enthropy_correction_pg=zad1_p_correction_pg(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb,om1_cychex,om1_eth);
enthropy_correction1=R.*log(p1_mix2./p1);
enthropy_without_correction=((integral(@(t) zad1_enthropy(t,tab1_moduls_cychex,M1_cychex),tem1_mix1,tem1_mix2)).*xa+(integral(@(t) zad1_enthropy(t,tab1_moduls_eth,M1_eth),tem1_mix1,tem1_mix2)).*xb)/1000-enthropy_correction1;
enthropy_correction2=R.*(xa.*log(xa)+xb.*log(xb));

enthropy_correction_pg_mixrule=zad1_p_correction_pg_enthropy_mixrule(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb,om1_cychex,om1_eth);
enthropy_correction_pg=zad1_p_correction_pg_enthropy(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb,om1_cychex,om1_eth);
enthropy_correction_rk=zad1_p_correction_rk_enthropy(p1_mix2,(tem1_mix2),pc_cychex,tc_cychex,pc_eth,tc_eth,xa,xb);

enthropy_sum=enthropy_without_correction-enthropy_correction2;

disp('Zmiana entropii mieszaniny bez uwzglêdniania poprawki: '+string(enthropy_correction2)+' J/mol*K');
disp('Poprawka ciœnieniowa na entropie wynosi:');
disp('Penga Robinsona regu³a mieszania: '+string(enthropy_correction_pg_mixrule)+'J/(mol*K)');
disp('Penga Robinsona regu³a Kaya     : '+string(enthropy_correction_pg)+'J/(mol*K)');
disp('Redlicha Kwonga regu³a Kaya     : '+string(enthropy_correction_rk)+'J/(mol*K)');

disp('Produkcja entropii: '+string((enthropy_sum-enthropy_correction_pg_mixrule-q1_str./273.15).*m1));

%zad2

disp(' ');
disp('Zadanie 2');
figure(2);
hold all
title('Ciep³o w³aœciwe przy sta³ym ciœnieniu etanolu w zale¿noœci temperatury');
xlabel('Temperatura [K]');
zlabel('Ciep³o w³aœciwe przy sta³ym ciœnieniu etanolu [J/K]');
ylabel('U³amek molowy');

tem_matrix=tem2_p:0.1:tem2_k;

%ethanol
%https://app.knovel.com/web/view/itable/show.v/rcid:kpYCPDCECD/cid:kt009ZN3EB/viewerType:eptble/root_slug:yaws-critical-property/url_slug:table-44-heat-capacity?b-q=material_or_substance_name%3A%22nitrogen%22%20AND%20specific_heat_mf%3A%5B%20*%20TO%20*%20%5D&b-facet-selected=item_type_nospace%3Aigraph&b-group-by=true&b-dsQuery=%3Fmn%3Dnitrogen%3Aand%3Apn%3Dspecific%20heat&b-search-type=tech-reference&b-sort-on=default&b-toc-cid=kpYCPDCECD&b-toc-root-slug=yaws-critical-property&b-toc-url-slug=table-44-heat-capacity&b-toc-title=Yaws%27%20Critical%20Property%20Data%20for%20Chemical%20Engineers%20and%20Chemists&start=0&columns=1,2,3,6,4,5,12,13,15,16,17,14,7,8,9,10,11&q=64-17-5
heat_constansts_ethanol=[238.308038502899
-2.38064232646184
1.33170564597485E-02
-3.19962372614431E-05
3.15051120535816E-08];

cp_eth=zad2_heat_eq1(tem_matrix,heat_constansts_ethanol,M2_ethanole);

%Water
%https://app.knovel.com/web/view/itable/show.v/rcid:kpYCPDCECD/cid:kt009ZN3DI/viewerType:eptble/root_slug:yaws-critical-property/url_slug:table-43-heat-capacity?b-q=material_or_substance_name%3A%22nitrogen%22%20AND%20specific_heat_mf%3A%5B%20*%20TO%20*%20%5D&b-facet-selected=item_type_nospace%3Aigraph&b-group-by=true&b-dsQuery=%3Fmn%3Dnitrogen%3Aand%3Apn%3Dspecific%20heat&b-search-type=tech-reference&b-sort-on=default&b-toc-cid=kpYCPDCECD&b-toc-root-slug=yaws-critical-property&b-toc-url-slug=table-43-heat-capacity&b-toc-title=Yaws%27%20Critical%20Property%20Data%20for%20Chemical%20Engineers%20and%20Chemists&start=0&columns=1,2,3,6,4,5,11,12,14,15,16,13,7,8,9,10&q=water
heat_constansts_water=[-22.4170167713592
0.876972155520207
-2.57039287293562E-03
2.48383325796128E-06
0];
cp_water=zad2_heat_eq1(tem_matrix,heat_constansts_water,M2_water);

%Mixture
%Table
surf(tem2_itc,xmol,matrix_from_itc);

interpzad2 = interp2(tem2_itc,xmol,matrix_from_itc,tem_matrix,x2_mol,'spline');
mix_correction=trapz(tem_matrix,interpzad2);

%clean enthalpy
heat_kg_mix=integral(@(t) zad2_heat_eq1(t,heat_constansts_ethanol,M2_ethanole),tem2_p,tem2_k).*x2_mol+integral(@(t) zad2_heat_eq1(t,heat_constansts_water,M2_water),tem2_p,tem2_k).*(1-x2_mol);
heat_kg_mix_with_corect_tab=mix_correction;

%total
total_heat_tab=heat_kg_mix_with_corect_tab.*m2;
disp('Ciep³o potrzebne do ogrzania strumienia: '+string(total_heat_tab/1000)+' kJ/s');

%zad3

disp(' ');
disp('Zadanie 3');
figure(3);
hold all;
title('Wykres logarytmu prê¿noœci pary nasyconej od temperatury');
xlabel('Temperatura [K]');
ylabel('lnP [mmHg]');

plot(tem3_zak,yaws3,'b');

%Drawing points
scatter(t3_vapour_temperature,log(p3_vapour_pressure),'xb');

%Computing constants
[a,b,c,error_anthonie_fit,res_matlab]=zad3_anthonie_find_constant(p3_vapour_pressure,t3_vapour_temperature,17.8417,-4796.1609,-40.4281);
disp('Sta³e wynosz¹:');
disp('A:  '+string(a));
disp('B: '+string(b));
disp('C: '+string(c));
disp('B³¹d przybli¿enia wynosi: '+string(error_anthonie_fit));

%Ploting fit curve and extarpolation for our range
plot(tem3_zak,(a+b./(tem3_zak+c)),'y');

legend('Date from Yaws','Date from Perry','Extrapolation');

%Pressure vapour
figure(6);
hold all
title('Wykres prê¿noœci pary nasyconej od temperatury');
xlabel('Temperatura [K]');
ylabel('P [Pa]');

semilogy([tem3,tem3],[0,Hg.*(exp(a+b./(tem3+c)))],'r');
scatter(tem3,Hg.*(exp(a+b./(tem3+c))),'or');
semilogy(tem3_zak,Hg.*exp(a+b./(tem3_zak+c)),'b');
semilogy([min(tem3_zak),tem3],[Hg.*(exp(a+b./(tem3+c))),Hg.*(exp(a+b./(tem3+c)))],'r');

vapour_pressure_Nma=Hg.*(exp(a+b./(tem3+c)));
disp('Prê¿noœæ pary nasyconej wynosi: '+string(vapour_pressure_Nma)+'Pa');

%Vapour heat
disp('Cip³o parowania wynosi: ')
%Carruth and Kobayashi
disp('Metoda Carruth-Kobayashi');
check_ck=tem3/tc3;
disp('Tempearatura zredukowana: '+string(check_ck));
if check_ck>0.6 && check_ck<1
    disp('Metoda Carruth-Kobayashi mo¿e zostaæ u¿yta');
    Vapour_heat_ck=zad3_vapour_heat_Carruth(tem3,tc3,acentric3);
disp(string(Vapour_heat_ck/1000)+' kJ/mol');
disp(string(Vapour_heat_ck/1000.*M3)+'  kJ/kg');
else
    disp('Metoda Carruth-Kobayashi NIE mo¿e zostaæ u¿yta');
end
%Clausiusa-Clapeyrona
disp('Metoda Clausiusa-Clapeyrona');
Vapour_heat_cc=zad3_vapour_heat_kc(tem3,b,c);
disp(string(Vapour_heat_cc/1000)+' kJ/mol');
disp(string(Vapour_heat_cc/1000.*M3)+'  kJ/kg');

l=zad3_redlich_kwong_for_cc(101315,tem3,pc3_for_cc,tc3_for_cc,a,b,c);

%Chart for dr Machniewski
figure(5);
hold all
plot(1./tem3_zak,(a+b./(tem3_zak+c)),'r');
scatter(1./t3_vapour_temperature,log(p3_vapour_pressure),'xb');
title('Wykres logarytmu prê¿noœci pary nasyconej od odwrotnoœci temperatury');
xlabel('Odwrotnoœæ temperatury [1/K]');
ylabel('lnP [mmHg]');

%Residua
[res3,sum_res3]=zad3_residua(p3_vapour_pressure,t3_vapour_temperature,a,b,c);
disp('Residua dla kolejnych punktów wynosz¹: ');
disp(res3');
disp('Suma residuów wynosi: '+string(sum_res3));

%zad4

disp(' ');
disp('Zadanie 4');
figure(7);
hold all

%http://app.knovel.com/web/view/itable/show.v/rcid:kpYCPDCECD/cid:kt00AAAI51/viewerType:eptble/root_slug:table-220-viscosity-of-aqueous-solution---inorganic-compounds-viscosity--a--bx--cx2--dx3--ex4/url_slug:table-220-viscosity-aqueous?q=material_or_substance_name%3A%227783-20-2%22%20AND%20dynamic_viscosity_mf%3A%5B%20*%20TO%20*%20%5D&b-q=material_or_substance_name%3A%227783-20-2%22%20AND%20dynamic_viscosity_mf%3A%5B%20*%20TO%20*%20%5D&sort_on=default&b-group-by=true&b-dsQuery=%3Fmn%3D7783-20-2%3Aand%3Apn%3Ddynamic%2520viscosity%26o%3De&dsQuery=%3Fmn%3D7783-20-2%3Aand%3Apn%3Ddynamic%2520viscosity%26o%3De&b-search-type=tech-reference&b-sort-on=default&scrollto=material_or_substance_name%3A%227783-20-2%22%20AND%20dynamic_viscosity_mf%3A%5B%20*%20TO%20*%20%5D&columns=1,2,3,4,5,6,7,11,8,9,12,13,14,10,15,16,17,18,19,10101010&start=0

c_p_computing=c4_salt.*M4/(c4_salt.*M4+1);
disp('Stê¿enie procentowe roztworu wynosi: '+string(c_p_computing*100)+'%');

%range [0.005 0.400]*100%
viscosity4_matrix_t20=[ 1.000238717864750
                        1.553360044484320
                       -0.298520311390234
                       15.567200889686400
                       -1.567200889686000];

viscosity4_t20=zad4_viscosity1(c_p_computing,viscosity4_matrix_t20);
disp('Lepkoœæ dla temperatury 20°C wynosi: '+string(viscosity4_t20)+'Pas');

%https://app.knovel.com/web/view/swf/show.v/rcid:kpICTNDPC4/cid:kt002VHUN2/viewerType:pdf/root_slug:international-critical?cid=kt002VHUN2&page=12&b-toc-cid=kpICTNDPC4&b-toc-root-slug=international-critical&b-toc-url-slug=viscosity-fluidity&b-toc-title=International%20Critical%20Tables%20of%20Numerical%20Data%2C%20Physics%2C%20Chemistry%20and%20Technology%20(1st%20Electronic%20Edition)
%kolumna t20 zostala uzupelniona z poprzedniego row0ia
tem3_range_matrix=[0,10,20,25,40,60,80]+K;
molality3_range_matrix=[0.25,0.5,1,2,2.5,4.5].*M4./([0.25,0.5,1,2,2.5,4.5].*M4+1);
Interntional_table_viscosity=[  NaN,    NaN,    1.0501, 1.053,  NaN,    NaN,    NaN;
                                NaN,    NaN,    1.092,  1.1,    1.112,  1.134,  1.162;
                                1.13,   1.18,   1.198,  1.209,  1.236,  1.269,  1.3;
                                NaN,    NaN,    1.452,  NaN,    1.52,   1.57,   1.61;
                                1.46,   1.55,   1.5999, 1.61,   1.67,   1.74,   NaN;
                                NaN,    NaN,    2.32,   NaN,    2.43,   2.5,    2.56];

%surf(tem3_range_matrix,molality3_range_matrix,Interntional_table_viscosity);

ylabel('Lepkoœæ roztworu');
title('Lepkoœæ roztwor w zalenoœci od temperatury');
xlabel('Temperatura [K]');
plot (tem3_range_matrix,Interntional_table_viscosity);
figure(8);
hold all
ylabel('Lepkoœæ roztworu');
title('Lepkoœæ roztwor w zalenoœci od molalnoœci');
xlabel('Molalnoœæ [mol/kg]');
plot (molality3_range_matrix,Interntional_table_viscosity);

%interpzad5 = interp2(tem3_range_matrix,molality3_range_matrix,Interntional_table_viscosity,tem4,c_p_computing,'spline')

%surf(tem3_range_matrix,molality3_range_matrix,interpzad5);

%zad5

disp(' ');
disp('Zadanie 5');

%Volume mol
volume_mol_ester=amount_c*element_c+amount_h*element_h+amount_o_double*element_o_duble+amount_o_ester_high*element_o_ester_high;
disp('Objêtoœæ molowa octanu n-butylu: '+string(volume_mol_ester)+'m^3/mol');

disp('Wspó³czynnik dyfuzyjnoœci obliczony:');

%Wilkego-Chenga 10% %Perry
diffusivity_temp_wc=zad5_Diffusivity_wc(tem5_p,associate5,M5,viscosity_tem5_p,volume_mol_ester);
diffusivity_temk_wc=zad5_Diffusivity_wc(tem5_k,associate5,M5,viscosity_tem5_k,volume_mol_ester);
disp('metod¹ Wilkego-Chenga:');
disp('T=10° : '+string(diffusivity_temp_wc)+' m^2/s');
disp('T=50° : '+string(diffusivity_temk_wc)+' m^2/s');

%Hayduka i Minhasa 14% https://ir.library.oregonstate.edu/xmlui/bitstream/handle/1957/33114/KittidachaWitoon2000.pdf?sequence=1
diffusivity_temp_hm=zad5_diffusivity_hm(volume_mol_ester,tem5_p,viscosity_tem5_p);
diffusivity_temk_hm=zad5_diffusivity_hm(volume_mol_ester,tem5_k,viscosity_tem5_k);
disp('metod¹ Hayduka i Minhasa:');
disp('T=10° : '+string(diffusivity_temp_hm)+' m^2/s');
disp('T=50° : '+string(diffusivity_temk_hm)+' m^2/s');

%Siddiqi i Lucasa 20%
diffusivity_temp_sl=zad5_diffusivity_sl(volume_mol_ester,tem5_p,viscosity_tem5_p);
diffusivity_temk_sl=zad5_diffusivity_sl(volume_mol_ester,tem5_k,viscosity_tem5_k);
disp('metod¹ Siddiqi i Lucasa:');
disp('T=10° : '+string(diffusivity_temp_sl)+' m^2/s');
disp('T=50° : '+string(diffusivity_temk_sl)+' m^2/s');















