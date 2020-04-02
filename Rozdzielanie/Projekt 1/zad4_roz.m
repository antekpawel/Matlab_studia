clear
clc

l=2.6;
t0=89+273;
ws=0.102;
ds=0.0026;
t1=18+273;


ct0=interp1([353 373],[148 180],t0);
ct1=interp1([283 293],[80 88],t1);

w=100*l/(100+ct0);
s0=ct0*l/(100+ct0);
s1=ct1*l/(100+ct0);

m=s0-s1;
mp=m/ws;
wpp=1+mp;

lam=fsolve(@(x) dzban(x,ds,wpp),0);

teta=0:0.5:5;
for i=1:max(size(teta))
    x1(i)=integral(@(te) (1+te.*lam./ds).^3.*4.*te.*exp(-2.*te),0,teta(i));
end
mian=integral(@(te) (1+te.*lam./ds).^3.*4.*te.*exp(-2.*te),0,Inf);
tab(:,3)=x1/mian;
tab(:,1)=teta;
tab(:,2)=ds+teta.*lam;
tab(:,4)=1-x1/mian;
figure(1)
clf
hold all
grid on
grid minor
hold all
plot(tab(:,2)*1000,tab(:,4))
ylabel('$$ 1-X_p [-] $$','Interpreter','latex')
xlabel('$$ d_s [mm] $$','Interpreter','latex')
function F=dzban(x,ds,wpp)

g=integral(@(te) (1+te.*x./ds).^3.*4.*te.*exp(-2.*te),0,Inf);

F=wpp-g;
end









