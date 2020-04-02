function wykresy_MW(t,rav,Rav)
clc

rpt=[10 20 40]*1e-9;

figure(1)
clf
hold all
plot(t*1e9,Rav(:),'r-')
% plot(t*1e9,Rav(:,2),'b-')
% plot(t*1e9,Rav(:,3),'g-')
grid on
grid minor
box on
xlabel('Czas t [ns]','interpreter','latex')
ylabel('\''Srednia szybko\''s\''c $R_{CO,avg}(t) \left[\frac{mol}{m^{2}\cdot s}\right] $','interpreter','latex')
legend({'$\overline{r}_{m}$ = 10 nm','$\overline{r}_{m}$ = 20 nm','$\overline{r}_{m}$ = 40 nm'},'interpreter','latex')

figure(2)
clf
hold all
plot(t*1e9,rav(:,1),'r-','linewidth',1.5)
plot(t*1e9,rav(:,2),'b-','linewidth',1.5)
plot(t*1e9,rav(:,3),'g-','linewidth',1.5)
grid on
grid minor
box on
xlabel('Czas t [ns]','interpreter','latex')
ylabel('Szybko\''s\''c konwersji $\overline{r}_{CO} \left[\frac{mol}{s}\right] $','interpreter','latex')
legend({'$\overline{r}_{m}$ = 10 nm','$\overline{r}_{m}$ = 20 nm','$\overline{r}_{m}$ = 40 nm'},'interpreter','latex')


s=10:1:40;
y=[Rav(end,1),Rav(end,2),Rav(end,3)];
z=rpt*1e9;

% f=@(x,a) x(1).*a.^2+x(2)*a+x(3);
f=@(x,a) x(1)*log(x(2)*a);


figure(3)
clf
% hold all
plot(rpt*1e9,[Rav(end,1),Rav(end,2),Rav(end,3)],'ko','MarkerFaceColor','r')
% plot(s,x(1)*s.^2+x(2)*s+x(3),'r--')
% plot(s,x(1)*log(x(2)*s),'r--')
grid on
grid minor
box on
xlabel('\''Sredni rozmiar mikro-por\''ow $ \overline{r}_{m}$ [nm]','interpreter','latex')
ylabel('\''Srednia szybko\''s\''c $R_{CO,avg} \left[\frac{mol}{m^{2}\cdot s}\right] $','interpreter','latex')
% axis([8 42 2.62e-3 2.82e-3])


y=[rav(end,1),rav(end,2),rav(end,3)];
Y=y;
X=z;

fun=@(x,a) x(1).*sqrt(x(2)*a);

x=lsqcurvefit(fun,[10 1],X,Y)




figure(4)
clf

% plot(s,x(1)*log(x(2)*s),'g--','LineWidth',1)
hold all
plot(rpt*1e9,[rav(end,1),rav(end,2),rav(end,3)],'ko','MarkerFaceColor','g')

grid on
grid minor
box on
xlabel('\''Sredni rozmiar mikro-por\''ow $ \overline{r}_{m}$ [nm]','interpreter','latex')
ylabel('Szybko\''s\''c konwersji CO $\overline{r}_{CO,avg} \left[\frac{mol}{s}\right] $','interpreter','latex')
% axis([8 42 9.8e-16 10.6e-16])


