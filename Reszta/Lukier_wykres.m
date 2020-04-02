


d=(66:125)/1e3;

for i=1:max(size(d))
    dp(i)=Lukier(14.61,50.58,0.842,d(i));
end

dp=dp'/1e5;

wykres(1)
plot(d,dp)


function dp=Lukier(tp,k,n,d)

%dane cieczy
ro=1500;
mk=0.556;


%sta³e
a=1./(2.*n+1);
b=2.*n./((n+1).*(2.*n+1));
c=2.*n.*n./((n+1).*(2.*n+1));

X=4.*tp./d;

%predkosc
v=4.*mk./pi./d.^2./ro;

%Kolejne cz³ony w równaniu
I  =4.*k./d;
II =(8.*v./d).^n;
III=((3.*n+1)./4./n).^n;
Suma=I.*II.*III;

%fun=@(press-Suma.*(1./(1-X./press)).*(1./(1-a.*X./press-b.*X.^2./press.^2-c.*X.^3./press.^3)).^n)
[dp]=fsolve(@(press) (press-Suma.*(1./(1-X./press)).*(1./(1-a.*X./press-b.*X.^2./press.^2-c.*X.^3./press.^3)).^n),1.5e5);
press=1e3:1e6;
dpp=press-Suma.*(1./(1-X./press)).*(1./(1-a.*X./press-b.*X.^2./press.^2-c.*X.^3./press.^3).^n);
plot(press,dpp)
axis([1 1e6 -1e-3 1e-3])
end