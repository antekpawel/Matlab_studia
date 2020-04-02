function [res,sum_res]=zad3_residua(ptab,ttab,a,b,c)

k=max(size(ptab));

%ptab and ttab is texperimental points

sum_res=0;
pom1=zeros(1,k);
pom2=zeros(1,k);
res=zeros(1,k);

ptab_m=log(ptab);

for i=1:k
    pom1(i)=ptab_m(i);
    pom2(i)=a+b./(ttab(i)+c);
    res(i)=abs(pom1(i)-pom2(i));
    sum_res=sum_res+abs(pom1(i)-pom2(i));
end