clc
clear

%% Dane
R = 8.314;
% Wiktor
k10 = 2.1833e+14;
E1 = 121000;
k20 = 1.7417e+26;
E2 = 182000;
Tmax = 361;
Tmin = 321;
CA0 = 1021;
beta = 1.31;
gamma = 0.021;
alfakon = 0.81;

%% Obliczenia
% Temperatura optymalna
alfa = [0:0.01:0.99];
Topt = (E2-E1)./(R.*log((k20.*(gamma+alfa).*E2)./(CA0.*k10.*(1-alfa).*(beta-alfa).*E1)));

figure(1)
hold all
grid on
grid minor
xlabel('{\alpha}_{A} [-]')
ylabel('T({\alpha}_{A}) [K]')
plot(alfa,Topt,'LineWidth',1)
plot([0 1],[Tmin Tmin],'--k','LineWidth',1)
plot([0 1],[Tmax Tmax],':k','LineWidth',2)
legend('T({\alpha}_{A})','T_{min}','T_{max}')

% Czas trwania procesu
alfamax = interp1(Topt,alfa,Tmax);
alfamin = interp1(Topt,alfa,Tmin);

fun1 = @(x) 1./((CA0.*k10.*exp(-E1./(R.*Tmax)).*(1-x).*(beta-x))-(k20.*exp(-E2./(R.*Tmax)).*(gamma+x)));
t1 = integral(fun1,0,alfamax);

% Robocze T(alfa)
% ((E2-E1)./(R.*log((k20.*(gamma+x).*E2)./(CA0.*k10.*(1-x).*(beta-x).*E1))))

fun2 = @(x) 1./((CA0.*k10.*exp(-E1./(R.*((E2-E1)./(R.*log((k20.*(gamma+x).*E2)./(CA0.*k10.*(1-x).*(beta-x).*E1)))))).*(1-x).*(beta-x))-(k20.*exp(-E2./(R.*((E2-E1)./(R.*log((k20.*(gamma+x).*E2)./(CA0.*k10.*(1-x).*(beta-x).*E1)))))).*(gamma+x)));
t2 = integral(fun2,alfamax,alfamin);

fun3 = @(x) 1./((CA0.*k10.*exp(-E1./(R.*Tmin)).*(1-x).*(beta-x))-(k20.*exp(-E2./(R.*Tmin)).*(gamma+x)));
t3 = integral(fun3,alfamin,alfakon);

tcalk = t1+t2+t3;

% t(alfa)

alfa01 = [0:0.0001:alfamax];
[m,n] = size(alfa01);
for i=1:n
t01(i) = integral(fun1,0,alfa01(i));
end

alfa02 = [alfamax:0.0001:alfamin];
[m,n] = size(alfa02);
for i=1:n
t02(i) = integral(fun2,alfamax,alfa02(i));
T02(i) = interp1(alfa,Topt,alfa02(i));
end

alfa03 = [alfamin:0.0001:alfakon];
[m,n] = size(alfa03);
for i=1:n
t03(i) = integral(fun3,alfamin,alfa03(i));
end

% Profil temperatury

figure(2)
hold all
grid on
grid minor
xlabel('t [s]')
ylabel('T [K]')
plot([0 t1],[Tmax Tmax],t1+t02,T02,[t1+t2 tcalk],[Tmin Tmin],'LineWidth',1)

% alfa(t)

figure(3)
hold all
grid on
grid minor
xlabel('t [s]')
ylabel('{\alpha}_{A} [-]')
plot(t01,alfa01,t1+t02,alfa02,t1+t2+t03,alfa03,'LineWidth',1)

disp('Poœrednie czasy i stopnie przereagowania, kiedy przestaj¹ dzia³aæ ograniczenia:')
disp('Czas t1, po którym zostanie osi¹gniête alfamax')
disp('Czas t2, po którym zostanie osi¹gniête alfamin')
disp('Czas t3, po którym zostanie osi¹gniête alfakon (docelowy stopieñ przemiany)')
disp('tcalk - ca³kowity czas procesu')