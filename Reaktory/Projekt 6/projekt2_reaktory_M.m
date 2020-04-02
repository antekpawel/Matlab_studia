R=1.905;                                     %[m] Promieñ reaktora
L=9.144;                                     %[m] D³ugoœæ reaktora
Cp=0.8374;                                    %[J/kgK] Ciep³o w³aœciwe
rho=1.6e-3;                                        %[kg/m3] Gêstoœæ mieszaniny
u=91440;                                %[m/s] Prêdkoœæ pozorna
rhob=1.443;                                   %[kg/m3] Gêstoœæ 
H=-69714;                                     %[J/kg]
kr=0.9353;                             %[J/s-m-K]
DR=4645;                              %[m2/s]
cA0=5e-5;                                     %[kg/m3]

Nm=8;
Nl=20;

dr=R/(Nm-1);
dz=L/(Nl-1);

f=zeros(Nm,Nl);
T=zeros(Nm,Nl);

f(:,1)=0;
T(:,1)=316.2+273.15;
T(Nm,:)=291+273.15;


% Warunki brzegowe

for l=1:Nl-1
    f(1,l+1)=f(1,l)+4*DR*dz/(u*(dr^2))*(f(2,l)-f(1,l))+ra(f(1,l),T(1,l))*rhob*dz/u;
    T(1,l+1)=T(1,l)+4*kr*dz/(rho*Cp*u*(dr^2))*(T(2,l)-T(1,l))-ra(f(1,l),T(1,l))*rhob*cA0*H*dz/(rho*Cp*u);
    for m=2:Nm-1
        f(m,l+1)=f(m,l)+DR*dz/u/dr^2*(1/(m-1)*((f(m+1,l)-f(m,l)))+f(m+1,l)-2*f(m,l)+f(m-1,l))+ra(f(m,l),T(m,l))*rhob*dz/u;
        T(m,l+1)=T(m,l)+kr*dz/(rho*Cp*u*(dr^2))*(1/(m-1)*(T(m+1,l)-T(m,l))+T(m+1,l)-2*T(m,l)+T(m-1,l))-ra(f(m,l),T(m,l))*rhob*cA0*H*dz/(rho*Cp*u);
    end
    f(Nm,l+1)=4/3*f(Nm-1,l+1)-1/3*f(Nm-2,l+1);
end
T=T';
f=f';
Ca=cA0*(1-f);
[x,y]=meshgrid(linspace(0,R,Nm),linspace(0,L,Nl));
size(T)
size(x)
size(y)

figure(1)
clf
hold all
surf(x,y,T-273.15,'edgecolor','none','facecolor','interp')
colormap(jet)
colorbar
xlabel('Promieñ R [m]')
ylabel('D³ugoœæ L [m]')

figure(2)
clf
hold all
surf(x,y,f,'edgecolor','none','facecolor','interp')
colormap(jet)
colorbar
xlabel('Promieñ R [m]')
ylabel('D³ugoœæ L [m]')

figure(3)
clf
hold all
surf(x,y,cA0*(1-f),'edgecolor','none','facecolor','interp')
colormap(jet)
colorbar
xlabel('Promieñ R [m]')
ylabel('D³ugoœæ L [m]')

function r=ra(f,T)

r=1223040*exp(-5389/T)*(1-f);
end

