function mi=zad4_viscosity1(cp,tab)
pom=tab(1)+tab(2).*cp+tab(3).*cp.*cp+tab(4).*cp.*cp.*cp+tab(5).*cp.*cp.*cp.*cp;
mi=pom./1000;