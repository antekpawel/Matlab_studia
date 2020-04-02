% function Intens_Proj1

e=0.3;
ep=0.51;
DL=1.2e-5;
KA=144.1;
ka=260;
DM=4.8e-5;
De=1.48e-5;
k=1.23e-4;

x=1;
u=0.1;
td=1;
dp=1e-3:1e-3:100e-3;

Bi=k.*dp./2./De;
K=1./(1./ka+1./k./KA);
fi=dp./2.*sqrt(ep.*K./De);
A0=fi.*coth(fi)-1;
G=12.*(1-e).*De.*Bi./e./dp.^2.*(1-Bi./(A0+Bi));
G0=sqrt(1+4.*DL.*G./u.^2);

A1=1./fi.*coth(fi)-(csch(fi)).^2;
Kpr=1+ka.^2.*KA./(ka+k.*KA).^2;
Gpr=1+3.*(1-e).*ep./2./e.*Bi.^2.*A1.*Kpr./(A0+Bi).^2;
Kbis=-2.*(ka.^2.*KA).^2./(ka+k.*KA).^3;
A2=1./fi.*(coth(fi)./fi.^2+csch(fi).^2./fi-2.*csch(fi).^2.*coth(fi));
Gbis=3.*(1-e).*ep./2./e.*Bi.^2./(A0+Bi).^3.*(-ep.*dp.^2./8./De.*(2.*A1.^2+(A0+Bi).*A2).*Kpr.^2+(A0+Bi).*A1.*Kbis);

mi0=exp(-u./2.*DL.*(G0-1).*x);
mi1=td./2+x./u./G0.*(1+3.*(1-e).*ep./2./e.*Bi.^2.*A1.*Kpr./(A0+Bi).^2);
mi2=td.^2./12+x./u.*(2.*DL./u.^2./G0.^3.*Gpr.^2-1./G0.*Gbis);

figure(1)
clf
hold all
plot(dp,1-mi0,'r-')
% plot(x./u,mi2,'b-')
% plot(1./u,mi2,'g-')
grid on
grid minor
box on
X='$$ d_p [m] $$';
% X='$$ \tau=\frac{L}{u} [s] $$' ;
% X='$$ \tau=\frac{L}{u} [s] $$';
Y='$ \alpha= 1-\mu_0 $';
% Y='$\tau_r=\mu_1 $ [s]';
% Y='$\sigma=\mu_2 [s^2] $';
podpis(X,Y)
