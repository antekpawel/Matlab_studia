%Funkcja do obliczania gêstoœci wody na podstawie danych z NIST'a
function ro=Ro(t)
I=importdata('Dane_woda.mat');
T=I(:,1);
RO=I(:,3);
ro=interp1(T,RO,t,'PCHIP');
end

