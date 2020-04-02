% function Zadani_2
clc
clear
Fl=0.1;           %m2
alfal=80;      %deg
DPl=0.46;       %At
DPl=DPl*98066.5;
q1=36.54;       %l/min  
q1=q1*1e-3/60;
q2=46.74; 
q2=q2*1e-3/60;
n1=1.4;        %1/min
n1=n1/60;
n2=2.85;
n2=n2/60;

s=0.98;
n=2.5;
n=n/60;
q=15.45;
q=q/3600;
alfa=100;
DP=0.506; %At
DP=DP*98066.5;

fun=@dzban;
x0=[0,0,0,0];
X=fsolve(fun,x0);
kl=X(1);
v0l=X(2);
v1=X(3);
v2=X(4);

antek=fsolve(@(u) powierzchnia(u,kl,DP,DPl,s,Fl,v0l,alfa,n,q),5)

function F=dzban(x)

Fl=0.1;           %m2
alfal=80;      %deg
DPl=0.46;       %At
DPl=DPl*98066.5;
q1=36.54;       %l/min  
q1=q1*1e-3/60;
q2=46.74; 
q2=q2*1e-3/60;
n1=1.4;        %1/min
n1=n1/60;
n2=2.85;
n2=n2/60;

F(1)=q1-n1*(x(3)-x(2));
F(2)=q2-n2*(x(4)-x(2));
F(3)=x(3)^2-x(2)^2-x(1)/n1*alfal/360;
F(4)=x(4)^2-x(2)^2-x(1)/n2*alfal/360;
%x1=kl
%x2=v0l
%x3=v1
%x4=v2
end
function G=powierzchnia(x,kl,DP,DPl,s,Fl,v0l,alfa,n,q)
    
    K=kl*((DP/DPl)^(1-s))*((x/Fl)^2);
    v0=v0l*x*((DP/DPl)^s)/Fl;
    v=(v0^2+K*alfa/(360*n))^0.5;
    qobl=n*(v-v0);
    G=abs(qobl-q);
end
    