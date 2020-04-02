function z=zad1_soave_2mix(p,t,pc1,tc1,pc2,tc2,x1,x2,om1,om2)

%sta³e
r=8.3144;u=2;w=-1;

om=x1.*om1+x2.*om2;

tc=x1.*tc1+x2.*tc2;

pc=x1.*pc1+x2.*pc2;

%Sk³adnik A
s=0.48+1.57*om-0.176*om^2;
a=0.42748*r^2*tc^2/pc*(1+s*(1-(t/tc)^0.5))^2;
b=0.08644*r*tc/pc;

%Obliczenie "duzych" wyrazów
A=a.*p./r^2./t^2;
B=b*p/r/t;

%pierwiastki rownania
z=roots([1,-1-B+u*B,A+w*B^2-u*B-u*B^2,-A*B-w*B^2-w*B^3]);