function mi=zad1_viscosity_wilke(mi1,mi2,y1,y2,m1,m2)
fi1=((1+((mi1./mi2).^0.5).*((m1./m2)^0.25)).^2)./(2.*sqrt(2).*(1+(m1./m2)^0.5));
fi2=((1+((mi2./mi1).^0.5).*((m2./m1)^0.25)).^2)./(2.*sqrt(2).*(1+(m2./m1)^0.5));
mi=(y1.*mi1)/(y1+y2.*fi2)+(y2.*mi2)/(y1.*fi1+y2);
