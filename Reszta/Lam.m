%Funkcja do obliczania ciepla wlasciwego wody na podstawie danych z NIST'a
function lam=Lam(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
LAM=I(:,13);
lam=interp1(T,LAM,t,'PCHIP');
end