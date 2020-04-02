function mi=zad1_viscosity_gas_perry(t,tab)
mi=tab(1).*t^tab(2)./(1+tab(3)./t+tab(4)./t./t);