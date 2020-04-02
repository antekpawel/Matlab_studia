%% Dane do zadania

R=0.10;
delta=0.002;
ro=900;
n=9/60;
M=2e-3;
% R=0.29;
% delta=0.005;
% ro=1900;
% n=9/60;
% M=3e-3;
%% Obliczenia 

mi=M*delta/pi^2/n/R^4;
ni=mi/ro;

save('zad3_mi.txt','mi');
save('zad3_ni.txt','ni');
disp(mi);
disp(ni);







































