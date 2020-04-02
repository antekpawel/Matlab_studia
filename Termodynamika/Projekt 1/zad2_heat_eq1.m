function cp=zad2_heat_eq1(t,tab,M)
cp=tab(1)+tab(2).*t+tab(3).*t.*t+tab(4).*t.*t.*t+tab(5).*t.*t.*t.*t;
cp=cp./M;

