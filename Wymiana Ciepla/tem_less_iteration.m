function [time]=tem_less_iteration(u1,u2)
time=fsolve(@(t) rownanie3(u1,u2,t),3225);

function te=rownanie1(u,t)
te=0;
a=4.444e-7;
X=0.09;


for i=1:19
    te=te+2.*(sin(u(1,i))./(u(1,i)+sin(u(1,i)).*cos(u(1,i)))).*exp(-(u(1,i).^2).*(a.*t./((X./2).^2)));
end

function te=rownanie2(u,t)
te=0;
a=4.444e-7;
X=0.09;

for i=1:19
    te=te+2.*((sin(u(1,i)))./(u(1,i)+sin(u(1,i)).*cos(u(1,i)))).*exp(-(u(1,i).^2).*(a.*t./((7*X./2).^2)));
end

function x=rownanie3(u1,u2,t)
teta=0.68164;

x=teta-(rownanie1(u1,t).^2).*rownanie2(u2,t);
