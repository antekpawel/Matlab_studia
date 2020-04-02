function ds=zad1_p_correction_pg_enthropy_mixrule(p,t,pc1,tc1,pc2,tc2,x1,x2,om1,om2)

%Coefficient
r=8.3144;u=2;w=-1;
%A
s1=0.37464+1.54226*x1.*om1-0.26992*x1.*om1.^2;
a1=0.45724*r^2*tc1^2/pc1*(1+s1*(1-(t/tc1)^0.5))^2;
b1=0.0778*r*tc1/pc1;


%B
s2=0.37464+1.54226*x1.*om2-0.26992*x1.*om2.^2;
a2=0.45724*r^2*tc2^2/pc2*(1+s2*(1-(t/tc2)^0.5))^2;
b2=0.0778*r*tc2/pc2;


%mix
a=(x1.*sqrt(a1)+x2.*sqrt(a2)).^2;
b=x1.*b1+x2.*b2;
%alpha

%Obliczenie "duzych" wyrazów
A=a.*p./r^2./t^2;
B=b*p/r/t;

%derrivetive
syms x;
y=(x1.*sqrt((0.45724*r^2*tc1^2/pc1*(1+s1*(1-(x/tc1)^0.5))^2))+x2.*sqrt((0.45724*r^2*tc2^2/pc2*(1+s2*(1-(x/tc2)^0.5))^2))).^2;
dy=diff(y) ;
D=subs(dy,t);
D=double(D);

%pierwiastki rownania
z=roots([1,-1-B+u*B,A+w*B^2-u*B-u*B^2,-A*B-w*B^2-w*B^3]);

%Enthalpy
ds=-r.*log(1-B./z(3))-r.*log(z(3))+D./b.*log((z(3)+B.*(1-sqrt(2)))./(z(3)+B.*(1+sqrt(2))));
