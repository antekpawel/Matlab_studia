function mi=zad1_viscosity_herning_zipperer(xa,xb,M1,M2,tc1,tc2,mi1,mi2)

%top
top=mi1.*xa.*sqrt(M1.*tc1)+mi2.*xb.*sqrt(M2.*tc2);

%bot
bot=xa.*sqrt(M1.*tc1)+xb.*sqrt(M2.*tc2);

%lepkoœæ
mi=top./bot;
