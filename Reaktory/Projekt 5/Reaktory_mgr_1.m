%% Pawe³ Antkowiak
%Reaktory mgr. Projekt 1

clear
close all
clc

R=8.314;
K=273.15;

%% Dane projetowe

u=2;
q=[1,3,10];
tau=0.5;
eta=2;

%logspace

%% Obliczenia pierwsze

l=tau*u;

z=0:0.01:l;

is=exp(-z./u./tau);
wykres(1)
plot(z,is,'b');
plot([max(z) max(z)],[0 1],'k--');
podpis('Odleglosc od wlotu [m]','Stopien segregacji [-]');
tabela1=logspace(-1,log10(l),10);
is2=exp(-tabela1./u./tau);
plot(tabela1,is2,'bx');


%%
for i=1:max(size(z))
    xa1(i)=fsolve(@(x) x^2+x.*(eta-1)-eta.*is(i),1);
end
for i=1:max(size(tabela1))
    fa(i)=1-fsolve(@(x) x^2+x.*(eta-1)-eta.*is2(i),1);
end

wykres(3)
plot(z,1-xa1,'b');
plot([max(z) max(z)],[0 1],'k--');
plot(tabela1,fa,'bx');
podpis('Odleglosc od wlotu [m]','Stopien przereagowania [-]');

%% Obliczenia drugie

xva=q./(q+1);
beta=q.*eta;


fs=beta./(1+beta);

for i=1:3
v(i,:)=xva(i).*((1./is)-1);
w(i,:)=(1-xva(i)).*(1./is-1);
end

for i=1:3
    for j=1:max(size(z))
        l1=1./v(i,j).*(0.001).^v(i,j);
        l2=1./w(i,j).*(0.001).^w(i,j);
        l3=integral(@(x) x.^(v(i,j)-1).*(1-x).^(w(i,j)-1),0.001,0.999);
        B(i,j)=l1+l2+l3;
    end
end
%B=B';
for i=1:3
    for j=1:max(size(z))
        Xaa(i,j)=1./xva(i).*integral(@(f) (f.*(1+beta(i))-beta(i)).*f.^(v(i,j)-1).*(1-f).^(w(i,j)-1)./B(i,j) ,fs(i),1-1e-10);
    end
end

wykres(4)
plot(z,Xaa)

%% kurwa
tab=[0.01000	0.00559	0.00261	0.00200
0.10000	0.05400	0.02600	0.01500
0.20000	0.10600	0.05100	0.03000
0.30000	0.15400	0.07700	0.04600
0.40000	0.20000	0.10100	0.06100
0.50000	0.24300	0.12600	0.07700
0.60000	0.28400	0.15000	0.09200
0.70000	0.32200	0.17400	0.10800
0.80000	0.35900	0.19800	0.12400
0.90000	0.39300	0.22100	0.14000
1.00000	0.42600	0.24400	0.15600];

wykres(5)
plot(tab(:,1),tab(:,2),'b');
plot([max(z) max(z)],[0 1],'k--');
plot(tab(:,1),tab(:,2),'bx');

plot(tab(:,1),tab(:,3),'r');
plot([max(z) max(z)],[0 1],'k--');
plot(tab(:,1),tab(:,3),'rx');

plot(tab(:,1),tab(:,4),'g');
plot([max(z) max(z)],[0 1],'k--');
plot(tab(:,1),tab(:,4),'gx');

podpis('Odleglosc od wlotu [m]','Stopien przereagowania [-]');













