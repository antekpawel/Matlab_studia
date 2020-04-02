function s=zad1_enthropy(t,tab,M)

pom=(tab(1)+tab(2).*t+tab(3).*t.*t+tab(4).*t.*t.*t+tab(5).*t.*t.*t.*t+tab(6).*t.*t.*t.*t.*t+tab(7).*t.*t.*t.*t.*t.*t)./t;

%Convert unit
s=pom*1000./M;