function l=zad3_vapour_heat_kc(t,b,c)

R=8.3144;

l=-(b.*R.*t.*t)./(c+t).^2;