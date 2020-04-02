function ds=zad1_p_correction_pg_enthropy(p,t,pc1,tc1,pc2,tc2,x1,x2,om1,om2)
%sta³e
r=8.3144;u=2;w=-1;

check_kay_rule=[tc1/tc2,pc1/pc2];
if check_kay_rule(1)<=2 && check_kay_rule(1)>=0.5 && check_kay_rule(2)<=2 && check_kay_rule(2)>=0.5

%Kay rule

s=0.37464+1.54226*(x1.*om1+x2.*om2)-0.26992*(x1.*om1+x2.*om2).^2;

tc=x1.*tc1+x2.*tc2;

pc=x1.*pc1+x2.*pc2;

%Coefficient

a=0.45724*r^2*tc^2/pc*(1+s*(1-(t/tc)^0.5))^2;
b=0.0778*r*tc/pc;

%alpha
alfa=1+s.*(1-(sqrt(t/tc))).^2;

%Obliczenie "duzych" wyrazów
A=a.*p./r^2./t^2.*alfa;
B=b*p/r/t;

%derrivetive
syms x ;
y=0.45724*r^2*tc^2/pc*(1+s*(1-(x/tc)^0.5))^2;
dy=diff(y) ;
D=subs(dy,t);
D=double(D);

%pierwiastki rownania
z=roots([1,-1-B+u*B,A+w*B^2-u*B-u*B^2,-A*B-w*B^2-w*B^3]);


%Enthalpy
ds=-r.*log(1-B./z(3))-r.*log(z(3))+D./b.*log((z(3)+B.*(1-sqrt(2)))./(z(3)+B.*(1+sqrt(2))));

else 
    disp('Nie mo¿na zastosowaæ regu³y Kaya :(');
    ds=0;
end