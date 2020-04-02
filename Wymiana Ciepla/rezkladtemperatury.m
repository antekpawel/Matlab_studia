function rezkladtemperatury(mi1,a1,l1,mi2,a2,l2,a,t,T0,Tp)

alin=l2/l1;
l1uz=l1.*sqrt(2);
l2uzy=sqrt(l2^2+l1^2);
k=max(size(mi1));
p=1;
for i=0:0.01:l2uzy/2
    sum1(p)=0;
    sum2(p)=0;
    for j=1:k
        sum1(p)=sum1(p)+(a1(j).*cos(mi1(j).*(2.*i./alin./l1uz)).*exp(-mi1(j)^2.*(a.*t.*4./l1^2)));
        sum2(p)=sum2(p)+(a2(j).*cos(mi2(j).*(2.*i./l2uzy)).*exp(-mi2(j)^2.*(a.*t.*4./l1^2)));
    end
    teta(p)=sum1(p)*sum1(p)*sum2(p);
    T(p)=teta(p).*(T0-Tp)+Tp;
    p=p+1;
end
figure(2)
hold all
plot(0:0.01:l2uzy/2,T,'r');
plot(0:-0.01:-l2uzy/2,T,'r');
grid on
xlabel('Odleg³oœæ [m]')
ylabel('Temperatura [K]')
title('Rozk³ad temperatury wzd³u¿ przek¹tnej')
%disp([0:0.01:l2uzy/2,T]');