%Funkcja do obliczania lepkoœci dynamicznej wody na podstawie danych z NIST'a
function mi=Mi(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
RO=I(:,12);
mi=interp1(T,RO,t,'PCHIP');
end

