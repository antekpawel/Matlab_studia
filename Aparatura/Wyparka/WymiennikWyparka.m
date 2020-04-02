%% Wymiennik 2 Wyparka
%Pawe³ Antkowiak
%Projekt aparatura

clear
close all
clc

g=9.81;% przyspieszenie ziemskie
K=273.15;
%% Dane

dw_pipe=[25 28 32 36 42].*1e-3;%[m]
delta_pipe=[2 2 2.5 3 3]*1e-3;%[m]

wew=1;%wybierz numer rury wew i zew
pipe_zew=800e-3;

t_wl_hot=100+K;
t_wyl_hot=36+K;
t_wl_cold=10+K;
t_wyl_cold=12+K;

m_czynnika=0.6944;
x_pary=0.3;
lams=58;
r=2.2567e6;

%% Obliczenia


dz_pipe=dw_pipe+2*delta_pipe;

pipe_wew=dz_pipe(wew);
pipe_wew_wew=dw_pipe(wew);

t_sr_hot=(t_wl_hot+t_wyl_hot)./2;
t_sr_cold=(t_wl_cold+t_wyl_cold)./2;

rho_hot=Ro(t_sr_hot);
rho_cold=Ro(t_sr_cold);

mi_hot=Mi(t_sr_hot);
mi_cold=Mi(t_sr_cold);

cp_hot=Cp(t_sr_hot);
cp_cold=Cp(t_sr_cold);

lam_hot=Lam(t_sr_hot);
lam_cold=Lam(t_sr_cold);

pr_hot=cp_hot.*mi_hot./lam_hot;
pr_cold=cp_cold.*mi_cold./lam_cold;

err=1;
err1=1;
err2=1;
N=316;
L2=10;
dT=1;

    Q=abs(m_czynnika.*integral(@(t) Cp(t),t_wl_hot,t_wyl_hot));
    Q1=m_czynnika.*x_pary.*r;
    m_wody=abs((Q+Q1)./integral(@(t) Cp(t),t_wl_cold,t_wyl_cold));

    v_w=4.*m_czynnika./pi./pipe_wew_wew.^2./rho_hot./N;
    v_z=m_wody./(pipe_zew.^2.*pi./4-pipe_wew.^2./4.*pi.*N)./rho_cold;

    Rew=pipe_wew_wew.*v_w.*rho_hot./mi_hot;
    dz=4.*(pipe_zew.^2./4.*pi-pipe_wew.^2./4.*pi.*N)./(pipe_zew.*pi+pipe_wew.*pi.*N);
    Rez=dz.*v_z.*rho_cold./mi_cold;
    
    Nuw=0.023.*(Rew^0.8).*(pr_hot^0.33);
    Nuz=0.027.*(Rez^0.8).*(pr_cold^0.33).*(mi_cold/Mi(18+K)).^0.14;
    
    alw=Nuw.*lam_hot./pipe_wew_wew;
    alz=Nuz.*lam_cold./dz;

    R=(1./pipe_wew_wew./pi./alw)+(log(pipe_wew./pipe_wew_wew)./2./pi./lams)+(1./pipe_wew./pi./alz);
    k=1/R;
    
while err2>1e-3
    while err1>1e-3
while err>1e-3
    
    alk=1.13.*(0.68204.^3.*917.01.^2.*g.*r./1.18246e-4.*dT./L2)^(1/4);

    R1=(1./pipe_wew_wew./pi./alk)+(log(pipe_wew./pipe_wew_wew)./2./pi./lams)+(1./pipe_wew./pi./alz);

    t_sr_log=((t_wl_hot-t_wyl_cold)-(t_wyl_hot-t_wl_cold))./log((t_wl_hot-t_wyl_cold)./(t_wyl_hot-t_wl_cold));
    L1=Q./N./k./t_sr_log;

    
    dT1=abs(Q1/alk/pi/pipe_wew_wew/L2);
    err=abs(dT1-dT);
    dT=dT1;
end
L21=(alk^4./1.13.^4./(0.68204.^3.*917.01.^2.*g.*r).*1.18246e-4.*dT).^(-1);
err1=abs(L21-L2);
L2=L21;
    end
        Nuw1=1.62*(Rew*pr_hot*pipe_wew_wew/L1).^0.33;
        err2=abs(Nuw1-Nuw);
        Nuw=Nuw1;
end
L=L1+L2;
%% Tworzenie tablic

dane_woda(1,:)=[t_sr_hot t_sr_cold];
dane_woda(2,:)=[rho_hot rho_cold];
dane_woda(3,:)=[mi_hot mi_cold];
dane_woda(4,:)=[cp_hot cp_cold];
dane_woda(5,:)=[lam_hot lam_cold];
dane_woda(6,:)=[pr_hot pr_cold];

%% Wyœwietlanie wyników

% disp('Œrednie wartoœci wody gor¹cej i zimnej')
% disp('Temperatura:')
% disp([t_sr_hot t_sr_cold])
% disp('Gêstoœæ:')
% disp([rho_hot rho_cold])
% disp('Lepkoœæ:')
% disp([mi_hot mi_cold])
% disp('Ciep³o w³aœciwe:')
% disp([cp_hot cp_cold])
% disp('Przewodnoœæ cieplna:')
% disp([lam_hot lam_cold])
% disp('Wspó³czynnik przewodzenia ciep³a:')
% disp([pr_hot pr_cold])

disp('Strumieñ ciep³a:')
disp(Q)
disp('Strumieñ masy czynnika chlodzacego:')
disp(m_wody)
disp('Predkoœæ czynnika wewnetrznego');
disp(v_w);
disp('Predkoœæ czynnika zewnetrznego');
disp(v_z);

disp('Raynolds czynnika wewnetrznego');
disp(Rew);
disp('Œrednica zastepcza');
disp(dz);
disp('Raynolds czynnika zewnetrznego');
disp(Rez);

disp('Nusselt czynnika wewnetrznego');
disp(Nuw);
disp('Nusselt czynnika zewnetrznego');
disp(Nuz);

disp('Alfa czynnika wewnetrznego');
disp(alw);
disp('Alfa czynnika zewnetrznego');
disp(alz);

disp('Wspolczynnik przenikania ciep³a');
disp(k);

disp('Srednia log');
disp(t_sr_log);

disp('Dlugosc wymiennika');
disp(L);

%% Zapis wyników
xlswrite('Tablica wartoœci wody',dane_woda);
disp(Mi(294.62))



%% Funkcje u¿ywane w rozwi¹zaniu problemu

%Funkcja do obliczania gêstoœci wody na podstawie danych z NIST'a
function ro=Ro(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
RO=I(:,3);
ro=interp1(T,RO,t,'PCHIP');
end

%Funkcja do obliczania lepkoœci dynamicznej wody na podstawie danych z NIST'a
function mi=Mi(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
RO=I(:,12);
mi=interp1(T,RO,t,'PCHIP');
end

%Funkcja do obliczania ciepla wlasciwego wody na podstawie danych z NIST'a
function cp=Cp(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
CP=I(:,9);
cp=interp1(T,CP,t,'PCHIP')*1000;
end

%Funkcja do obliczania ciepla wlasciwego wody na podstawie danych z NIST'a
function lam=Lam(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
LAM=I(:,13);
lam=interp1(T,LAM,t,'PCHIP');
end

