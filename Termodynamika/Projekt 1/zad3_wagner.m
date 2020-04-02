function p=zad3_wagner(t,tc,pc,a1,a2,a3)

%Reduct temperature
tr=t./tc;

%Eqation
pr=exp((a1.*(1-tr)+a2.*(1-tr)^(4./3)+a3.*(1-tr).*(1-tr).*(1-tr))./tr);

%Real pressure
p=pr*pc;