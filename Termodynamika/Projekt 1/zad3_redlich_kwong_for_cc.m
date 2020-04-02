function l=zad3_redlich_kwong_for_cc(p,t,pc,tc,a,b,c)

%t=temperatura [K], p=ciœnienie [Pa]

r=8.3144;u=1;w=2;

a1=0.42748*r^2*tc^2.5/pc/t^0.5;
b1=0.0866*r*tc/pc;

A=a1*p/r^2/t^2;
B=b1*p/r/t;

z=roots([1,-1-B+u*B,A+w*B^2-u*B-u*B^2,-A*B-w*B^2-w*B^3]);

syms x
y=exp(a-b./((1./x)-c));
dy=diff(y);
pom=subs(dy,t);
pom1=double(pom);

l=-r.*(z(3)-z(1)).*pom1;