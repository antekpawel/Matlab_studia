%Funkcja do obliczania ciepla wlasciwego wody na podstawie danych z NIST'a
function cp=Cp(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
CP=I(:,9);
cp=interp1(T,CP,t,'PCHIP')*1000;
end

