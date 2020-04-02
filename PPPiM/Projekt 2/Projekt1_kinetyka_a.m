%% Pawe³ Antkowiak
%Projekt 2 zad1
%Kinetyka

clear
clc
close all
%% Dane
a=3e-2;%[m]
b=5e-2;%[m]
c=9e-2;%[m]
wp=0.3;%[-]
Dab=1e-9;%[m2/s]
wr=0.05;
wk=0.15;
Bi=4;





%% Obliczenia
Cp=wp/(1+wp);
Ck=wk/(1+wk);
Cr=wr/(1+wr);

t=0:1e3:4e5;

mi=Mip(Bi,10);
ap=Anp(Bi,mi);


y1=0;
for i = 1:max(size(ap))
    y1=y1+ap(i).*exp(-mi(i).^2.*Dab.*t./(c/2).^2);
end

y2=0;
for i = 1:max(size(ap))
    y2=y2+ap(i).*exp(-mi(i).^2.*Dab.*t./(b/2).^2);
end
y=y1.*y2;

c1=y.*(Cp-Cr)+Cr;
w=c1./(1-c1);
hold all
plot(t,w);
chuj1=interp1(w,t,wk);
k=integral(@(x) dzban(chuj1,x),0,a/2)*2/a
disp(k);
y=k
c2=y.*(Cp-Cr)+Cr
w2=c2./(1-c2)















function y=dzban(t,x)
c=9e-2;%[m]
Dab=1e-9;%[m2/s]
Bi=50;
mi=Mip(Bi,10);
ap=Anp(Bi,mi);
y1=0;
for i = 1:max(size(ap))
    y1=y1+ap(i).*exp(-mi(i).^2.*Dab.*t./(c/2).^2).*cos(mi(i)*x*2/c);
end
y=y1;

end

function y=dzban1(t,x)
b=5e-2;%[m]
Dab=1e-9;%[m2/s]
Bi=50;
mi=Mip(Bi,10);
ap=Anp(Bi,mi);
y1=0;
for i = 1:max(size(ap))
    y1=y1+ap(i).*exp(-mi(i).^2.*Dab.*t./(b/2).^2).*cos(mi(i)*x*2/b);
end
y=y1;

end
