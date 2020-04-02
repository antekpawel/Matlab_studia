function z=zad1_soave_2mix_mixrule(p,t,pc1,tc1,pc2,tc2,x1,x2,om1,om2)

%sta³e
r=8.3144;u=2;w=-1;

%Sk³adnik A
s1=0.48+1.57*om1-0.176*om1^2;
a1=0.42748*r^2*tc1^2/pc1*(1+s1*(1-(t/tc1)^0.5))^2;
b1=0.08644*r*tc1/pc1;

%B
s2=0.48+1.57*om2-0.176*om2^2;
a2=0.42748*r^2*tc2^2/pc2*(1+s2*(1-(t/tc2)^0.5))^2;
b2=0.08644*r*tc2/pc2;

%
a=(x1.*sqrt(a1)+x2.*sqrt(a2)).^2;
b=x1.*b1+x2.*b2;

%Obliczenie "duzych" wyrazów
A=a.*p./r^2./t^2;
B=b*p/r/t;

%pierwiastki rownania
z=roots([1,-1-B+u*B,A+w*B^2-u*B-u*B^2,-A*B-w*B^2-w*B^3]);