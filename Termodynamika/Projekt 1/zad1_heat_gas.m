function cp=zad1_heat_gas(t,tab)

pom=tab(1)+tab(2).*t+tab(3).*t.*t+tab(4).*t.*t.*t+tab(5).*t.*t.*t.*t+tab(6).*t.*t.*t.*t.*t+tab(7).*t.*t.*t.*t.*t.*t;

%Convert unit
cp=pom*1000;