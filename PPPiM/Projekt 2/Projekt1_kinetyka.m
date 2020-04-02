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
wp=0.32;%[-]
Dab=1e-9;%[m2/s]
wr=0.05;
wk=0.27;
Bi=50;




%% Obliczenia
Cp=wp/(1+wp);
Ck=wk/(1+wk);
Cr=wr/(1+wr);

t=0:1e3:3e5;

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
%plot(t,y);





t=chuj1;
y21=0;
for i = 1:max(size(ap))
    y21=y21+ap(i).*exp(-mi(i).^2.*Dab.*t./(a).^2).*cos(mi(i)/2);
end
y21=1;
y22=0;
for i = 1:max(size(ap))
    y22=y22+ap(i).*exp(-mi(i).^2.*Dab.*t./(b).^2).*cos(mi(i)/2);
end

y23=0;
for i = 1:max(size(ap))
    y23=y23+ap(i).*exp(-mi(i).^2.*Dab.*t./(c).^2).*cos(mi(i)/2);
end
y23=1;
yk=y21.*y22.*y23;

c2=yk.*(Cp-Cr)+Cr;
w2=c2./(1-c2);
figure(2)
hold all
plot(t,w2);
plot(t,yk);





