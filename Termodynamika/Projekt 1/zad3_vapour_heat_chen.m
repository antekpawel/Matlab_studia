function l=zad3_vapour_heat_chen(t,tw,tc,pc)

%constanc
r=8.3144;bar1=100000;

%Reduced temperature
tbr=t./tw;

%Heat
l=r.*tc.*((3.978.*tbr-3.958+1.555.*log(pc./bar1))./(1.07-tbr));
