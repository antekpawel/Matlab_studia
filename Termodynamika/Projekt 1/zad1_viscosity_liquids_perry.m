function mi=zad1_viscosity_liquids_perry(t,tab)
mi=exp(tab(1)+tab(2)./t+tab(3).*log(t)+tab(4).*t^(tab(5)));