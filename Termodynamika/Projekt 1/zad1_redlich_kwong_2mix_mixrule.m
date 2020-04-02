function z=zad1_redlich_kwong_2mix_mixrule(p,t,pc1,tc1,pc2,tc2,x1,x2)
%t=temperatura [K], p=ciœnienie [Pa]
r=8.3144;u=1;w=0;

%
a1=0.42748*r^2*tc1^2.5/pc1/t^0.5;
b1=0.0866*r*tc1/pc1;

%
a2=0.42748*r^2*tc2^2.5/pc2/t^0.5;
b2=0.0866*r*tc2/pc2;

%
a=(x1.*sqrt(a1)+x2.*sqrt(a2)).^2;
b=x1.*b1+x2.*b2;

%
A=a.*p./r^2./t^2;
B=b.*p./r./t;

z=roots([1,-1-B+u.*B,A+w.*B^2-u.*B-u.*B.*B,-A.*B-w.*B.*B-w.*B.*B.*B]);