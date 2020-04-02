function ds=zad1_p_correction_rk_enthropy(p,t,pc1,tc1,pc2,tc2,x1,x2)
%t=temperatura [K], p=ciœnienie [Pa]
r=8.3144;u=1;w=0;

tc=x1.*tc1+x2.*tc2;

pc=x1.*pc1+x2.*pc2;

%
a=0.42748*r^2*tc^2.5/pc/t^0.5;
b=0.0866*r*tc/pc;

%
A=a.*p./r^2./t^2;
B=b.*p./r./t;

z=roots([1,-1-B+u.*B,A+w.*B^2-u.*B-u.*B.*B,-A.*B-w.*B.*B-w.*B.*B.*B]);

ds=r.*(A./B./2.*log(1+B./z(3))-log(z(3)-B./z(3)));