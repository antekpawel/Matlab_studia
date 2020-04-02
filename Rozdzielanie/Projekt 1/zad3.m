% function zad3
clear
clc
M=19.9;
t0=89+273;
Ws=1.02;
t1=18+273;

ass=0.5:0.5:4.5;
ds=[0 0.038 0.146 0.309 0.5 0.691 0.854 0.942 1];
WP=1-ds;
c0=53.713;
c1=33.41;
% c0=54.849;
% c1=31.311;
W=100*M/(100+c0);
S0=c0*M/(100+c0);

S1=c1*M/(100+c0);
M1=S1+W;

m=S0-S1;
mp=m/1;
Wpp=1+mp;

dd=fsolve(@(x) dzban3(x,Wpp,WP,ass),0)



figure(1)
clf
hold all
grid on
grid minor
hold all
y=0:1e-3:1;

plot(interp1(ds,ass,y,'PCHIP'),y,'b')
plot(interp1(ds,ass+dd,y,'PCHIP'),y,'r')

ylabel('$$ X_p $$','Interpreter','latex')
xlabel('$$ d_s $$','Interpreter','latex')
legend('Szczepionka','Produkt','Location','Northwest')

DELTA=dd;

ass=1:0.5:5;
ds=[0 0.038 0.146 0.309 0.5 0.691 0.854 0.942 1];
    
dzban2=(1+DELTA./ass).^3;
pole(1)=-trapz(WP,dzban2);

ass=1:0.5:4.5;
ds=[0 0.038 0.146 0.309 0.5 0.691 0.854 0.942];
WP=1-ds;
dzban2=(1+DELTA./ass).^3;
pole(2)=-trapz(WP,dzban2);

ass=1:0.5:4;
ds=[0 0.038 0.146 0.309 0.5 0.691 0.854];
WP=1-ds;
dzban2=(1+DELTA./ass).^3;
pole(3)=-trapz(WP,dzban2);

ass=1:0.5:3.5;
ds=[0 0.038 0.146 0.309 0.5 0.691];
WP=1-ds;
dzban2=(1+DELTA./ass).^3;
pole(4)=-trapz(WP,dzban2);

ass=1:0.5:3;
ds=[0 0.038 0.146 0.309 0.5];
WP=1-ds;
dzban2=(1+DELTA./ass).^3;
pole(5)=-trapz(WP,dzban2);

ass=1:0.5:2.5;
ds=[0 0.038 0.146 0.309];
WP=1-ds;
dzban2=(1+DELTA./ass).^3;
pole(6)=-trapz(WP,dzban2);

ass=1:0.5:2;
ds=[0 0.038 0.146];
WP=1-ds;
dzban2=(1+DELTA./ass).^3;
pole(7)=-trapz(WP,dzban2);

ass=1:0.5:1.5;
ds=[0 0.038];
WP=1-ds;
dzban2=(1+DELTA./ass).^3;
pole(8)=-trapz(WP,dzban2);


pole(9)=0;

pom=pole(1);
dzban4=1-pole./pom

function F=dzban(x,ds,wpp,ass)

g=integral(@(u) (1+x./interp1(ds,ass,u,'PCHIP')).^3,0,1);

F=abs(g-wpp);

end

function F=dzban3(x,Wpp,WP,ass)
dzban2=(1+x./ass).^3;
pole=-trapz(WP,dzban2);
F=abs(pole-Wpp);

end


