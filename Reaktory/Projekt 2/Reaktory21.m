%% Pawe³ Antkowiak
clear
close all
clc

%% Dane projektowe

k0=9.4e5;                           % [dm3/mol*s2]
ea=17e3;                            % [J/mol]
h=-139e3;                           % [J/mol]
cp=4.19e3;                          % [J/kg*K]
ro=0.8;                             % [kg/dm3]
ca0=23;                             % [mol/dm3]
cb0=4.6;                            % [mol/dm3]
vb0=5.7;                            % [dm3]
va0=1.5;                            % [dm3]
T1=22+273.15;                       % [K]
t1=8.5;                             % [s]
t2=107;                             % [s]

R=8.314;                            % [J/mol*K]
%% Stêzenie alfa=0.999
a=0.999;
cbx=(cb0*vb0)/(va0+vb0)*(1-a);
%% Obliczenia

q1=va0/t1;
k2=k0*exp(-ea/(R*T1));

IC=[0 cb0];

[tp,cp]=ode45(@(t,c) dzban1(q1,ca0,c,k2,t1,vb0,va0,t),[0 t1],IC);

%% Funkcje

function F=vol(t,vb0,va0,ti)
F=vb0+va0.*t./ti;
end
function F=dzban1(q,ca0,c,k2,t1,vb0,va0,t)

F(1)=q./vol(t,vb0,va0,t1).*(ca0-c(1))-k2.*c(1).*c(2);
F(2)=-q./vol(t,vb0,va0,t1).*c(2)-k2.*c(1).*c(2);
end









