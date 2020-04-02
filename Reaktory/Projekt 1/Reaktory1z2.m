clc
clear
%% Dane

% %Wiktor
% t = [10 30 50 70 100];
% ca0 = 2.5;
% ca = [2.05 1.42 1.02 0.75 0.49];
% alfa = [0.178 0.430 0.592 0.701 0.804];
% n = [1.01:0.001:1.4];

%Antke
t = [10 30 50 70 100];
ca0 = 3;
ca = [2.04 1.21 0.84 0.64 0.46];
% alfa = [0.083 0.241 0.388 0.524 0.704];
n = 1.6:0.01:2.1;

% %Elza
% t = [10 30 50 70 100];
% ca0 = 3.5;
% ca = [3.09 2.36 1.76 1.26 0.71];
% alfa = [0.118 0.325 0.498 0.639 0.798];
% n = [0.1:0.001:0.99];
% 
% %Kasia
% t = [10 30 50 70 100];
% ca0 = 4;
% ca = [3.29 2.20 1.46 0.97 0.51];
% alfa = [0.179 0.449 0.634 0.758 0.873];
% n = [0.1:0.001:0.999 1.001:0.001:1.4];

%% Czarna magia

% fsolve(@(u) dzban(n,t,ca0,ca),1)
% 
% function F=dzban(n,t,ca0,ca)
% 
%     k=sum((ca^(1-n))-(ca0^(1-n)))./(t*(n-1));
% 
% sr=sum(k)/5;
% F=sum((k-sr).^2)/5;
% 
% 
% end
[M,N] = size(n);

for i=1:N
    for j=1:5
    k(i,j)=((ca(j).^(1-n(i)))-(ca0.^(1-n(i))))./(t(j).*(n(i)-1));
    end
end

ksr=sum(k,2)./5;

for i=1:N
    for j=1:5
        rozn(i,j)=(k(i,j)-ksr(i)).^2;
    end
end

vark=sum(rozn,2)./5;

% Wykres
figure
plot(n,vark)
grid on
grid minor
xlabel('n [-]')
ylabel('var k ')

% Szukanie n
z=1;
while vark(z+1,1) < vark(z,1)
    nmin=n(1,z);
    z=z+1;
end

% Obliczenie k
kszukane=ksr(z,1)