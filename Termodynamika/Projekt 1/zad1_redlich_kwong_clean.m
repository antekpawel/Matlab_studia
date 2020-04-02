function z=zad1_redlich_kwong_clean(p,t,pc,tc)

%t=temperatura [K], p=ciœnienie [Pa]

r=8.3144;u=1;w=2;

a=0.42748*r^2*tc^2.5/pc/t^0.5;
b=0.0866*r*tc/pc;

A=a*p/r^2/t^2;
B=b*p/r/t;

z=roots([1,-1-B+u*B,A+w*B^2-u*B-u*B^2,-A*B-w*B^2-w*B^3]);
