clear
clc

d=10e-3;
Re=1000;
mi=1.0016e-3;
ro=1e3;
R=d/2;


r=0:1e-4:R;
u=Re*mi/(ro*d);
Q=pi*(d^2)*0.25*u;
dpl=8*Q*mi/(pi*R^4);
vx=dpl*(R^2)/(4*mi).*(1-(r.^2)/(R^2));

figure(1);
grid on;
grid minor
hold all
plot(vx,r,'r');
plot([0 max(vx)],[-max(r) -max(r)],'b');
%legend('Œcianka rury','Profil prêdkoœci','Location','Best')


plot(vx,-r,'r');
plot([0 max(vx)],[max(r) max(r)],'b');
xlabel('Prêdkoœæ cieczy [m/s]');
ylabel('Promieñ rury [m]');
title('Rozk³ad prêdkoœci dla Re='+string(Re));
axis([0 max(vx) -8e-3 8e-3]);

plot([0 max(vx)],[0 0],'-.','Color',[0 102 0]/255);

saveas(gcf,'F:\Pulpit\Nauka\Studia\III rok\PPPiM\Projekt1\TemperatureMeasurement_SCS\lam','epsc')

tab2(:,1)=r;
tab2(:,2)=vx;





%% B
Re=8e5;
y=R-r;

u=Re*mi/(ro*d);
lam=0.0032+0.221/Re^0.237;%Blasius
dpl=lam*ro*u^2/d/2;
tw=d*dpl/4;
u0=sqrt(tw/ro);
ypl=y*u0*ro/mi;

for i=1:max(size(ypl))
    if ypl(i)<=5
        v(i)=u0*ypl(i);
    elseif ypl(i)<=30 && ypl(i)>5
        v(i)=u0*(5*log(ypl(i))-3.05);
    else
        v(i)=u0*(2.5*log(ypl(i))+5.5);
    end
end

figure(2)
grid on;
grid minor
xlabel('Prêdkoœæ cieczy [m/s]');
ylabel('Promieñ rury [m]');
title('Rozk³ad prêdkoœci dla Re='+string(Re));
hold all
axis([0 max(v) -8e-3 8e-3]);
text(1,-4e-3,'Obszar podwarstwy laminarnej','Rotation',90);
text(4,-4e-3,'Obszar warstwy przejœciowej','Rotation',90);
text(25,-4e-3,'Obszar rdzenia turbulentnego','Rotation',90);
plot([0 max(v)],[-max(y) -max(y)],'b');
plot(v,y-max(y),'r');
%legend('Œcianka rury','Profil prêdkoœci','Location','Best')

plot([0 max(v)],[0 0],'-.','Color',[0 102 0]/255);
plot([0 max(v)],[max(y) max(y)],'b');
plot(v,-y+max(y),'r');

vlam=interp1(ypl,v,5);
vprz=interp1(ypl,v,30);
vtur=max(v);

disp(vlam);
disp(vprz);
disp(vtur);

plot([vlam vlam],[-max(y) max(y)],'k--');
plot([vprz vprz],[-max(y) max(y)],'k--');


saveas(gcf,'F:\Pulpit\Nauka\Studia\III rok\PPPiM\Projekt1\TemperatureMeasurement_SCS\burz','epsc')
tab(:,1)=r;
tab(:,2)=ypl;
tab(:,3)=v;















