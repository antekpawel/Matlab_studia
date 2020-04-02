function [Rav rav t]=MW_proj_1

%% Dane projektowe
clc
global k1 Deff gamma em rm rpt R M

j=1;
jj=1;

dpm=0.5e-6;             % [m]
rm=[10 20 40]*1e-9;     % [m]
em=0.25;                % [-]
tm=2.2;                 % [-]
rpt=[5 10 20]*1e-9;     % [m]
fipt=0.05;              % [u³. obj]
Ts=600;                 % [K]
ys0=0.001;              % [u³ mol.]

R=8.314;                % [J/mol/K]
M=28e-3;                % [kg/mol]

tgrid=300;


Dk=2*rm(j)*sqrt(8*R*Ts/pi/M)/3;
Deff=em*Dk/tm;
k1=7.673e10*exp(-12556/Ts);
Vm=pi*dpm^3/6;
gamma=3*fipt./rpt(jj);
A_Pt=Vm.*gamma;


r=linspace(0,dpm/2,20);
t=linspace(0,5e-9,tgrid);

m = 2;
sol = pdepe(m,@pdefun,@IC,@BC,r,t);

u(:,:,j) = sol(:,:,1);
RV(:,:,j)=k1*u(:,:,j)*gamma;
jj=2;
j=1;
Dk=2*rm(j)*sqrt(8*R*Ts/pi/M)/3;
Deff=em*Dk/tm;
k1=7.673e10*exp(-12556/Ts);
Vm=pi*dpm^3/6;
gamma=3*fipt./rpt(jj);
A_Pt=Vm.*gamma;


r=linspace(0,dpm/2,20);
t=linspace(0,5e-9,tgrid);

m = 2;
sol = pdepe(m,@pdefun,@IC,@BC,r,t);

u(:,:,jj) = sol(:,:,1);
RV(:,:,jj)=k1*u(:,:,jj)*gamma;
jj=3;
j=1;
Dk=2*rm(j)*sqrt(8*R*Ts/pi/M)/3;
Deff=em*Dk/tm;
k1=7.673e10*exp(-12556/Ts);
Vm=pi*dpm^3/6;
gamma=3*fipt./rpt(jj);
A_Pt=Vm.*gamma;


r=linspace(0,dpm/2,20);
t=linspace(0,5e-9,tgrid);

m = 2;
sol = pdepe(m,@pdefun,@IC,@BC,r,t);

u(:,:,jj) = sol(:,:,1);
RV(:,:,jj)=k1*u(:,:,jj).*gamma;
% 
% % nazwa='Profil_rm'+string(rm(j)*1e9)'+'_rpt'+string(rpt(jj)*1e9)+'.avi';
% % nazwa=convertStringsToChars(nazwa);
% % 
% % writerObj=VideoWriter(nazwa);
% % writerObj.Quality = 100;
% % open(writerObj)
% % 
% % figure(1)
% % clf
% for i=1:tgrid
%         hold all
%         plot(r.*1e9,u(i,:,1),'r-','LineWidth',1.5)
%         plot(r.*1e9,u(i,:,2),'b-','LineWidth',1.5)
%         plot(r.*1e9,u(i,:,3),'g-','LineWidth',1.5)
% %         plot(r.*1e9,u(i,:,4),'y-','LineWidth',1.5)
%         xlabel('Promie\''n r [nm]','interpreter','latex')
%         ylabel('St\c{e}\.zenie $y_{CO}(r,t)\,[-] $','interpreter','latex')
%         grid on
%         grid minor
%         box on
%         axis([0 250 0.0024 ys0])
%         text(50,0.00295,'Czas t = '+string(t(i)*1e9)+' ns','interpreter','latex')
%         legend({'$\overline{r}_{m} =\,10\,nm,\,r_{Pt} =\,5\,nm$'...
%             '$\overline{r}_{m} =\,20\,nm,\,r_{Pt} =\,10\,nm$'...
%             '$\overline{r}_{m} =\,40\,nm,\,r_{Pt} =\,20\,nm$'},'interpreter','latex','Location','SouthEast')
%         
%         F = getframe(figure(1));
%         writeVideo(writerObj,F);
%         clf
% end
% close(writerObj)
% 
figure(2)
clf
hold all
plot(r.*1e9,RV(end,:,1),'r-','LineWidth',1.5)
% plot(r.*1e9,RV(end,:,2),'b-','LineWidth',1.5)
% plot(r.*1e9,RV(end,:,3),'g-','LineWidth',1.5)
xlabel('Promie\''n r [nm]','interpreter','latex')
ylabel('Szybko\''s\''c reakcji $R_{CO,V}(r)\,[\frac{mol}{m^3\cdot s}] $','interpreter','latex')
grid on

grid minor
box on
% % axis([0 250 0.0024*k1*gamma ys0*k1*gamma])
legend({'$\overline{r}_{m} =\,10\,nm,\,r_{Pt} =\,5\,nm$'...
'$\overline{r}_{m} =\,10\,nm,\,r_{Pt} =\,10\,nm$'...
'$\overline{r}_{m} =\,10\,nm,\,r_{Pt} =\,20\,nm$'},'interpreter','latex','Location','SouthEast')
     


% nazwa='Szybkosc_rm'+string(rm(j)*1e9)+'_rpt'+string(rpt(jj)*1e9)+'.avi';
% nazwa=convertStringsToChars(nazwa);
% writerObj=VideoWriter(nazwa);
% writerObj.Quality = 100;
% open(writerObj)

% figure(2)
% for i=1:tgrid
%         plot(r.*1e9,k1*u(i,:,j),'b-','LineWidth',1.5)
%         xlabel('Promie\''n r [nm]','interpreter','latex')
%         ylabel('Szybko\''s\''c $R_{CO}(r,t) \left[\frac{mol}{m^{2}\cdot s}\right] $','interpreter','latex')
%         grid on
%         grid minor
%         box on
%         axis([0 250 0.0025 k1*ys0]) 
%         text(150,0.06625,'Czas t = '+string(t(i)*1e9)+' ns','interpreter','latex')
%         text(50,0.06935,'$\overline{r}_{m}$ = '+string(rm(j)*1e9)+' nm','interpreter','latex')
%         text(50,0.06915,'$r_{Pt}$ = '+string(rpt(jj)*1e9)+' nm','interpreter','latex')
%         
%         F = getframe(figure(2));
%         writeVideo(writerObj,F);
%         
% end
% close(writerObj)
% 

% 
figure(3)
clf
hold all
% surf(r*1e9,t*1e9,u(:,:,3),'edgecolor','none','facecolor','interp')
xlabel('Promieñ r [nm]')
ylabel('Czas t [ns]')
zlabel('St\c{e}\.zenie $y_{CO}(r,t) [-] $')
grid on
colormap(hot)
colorbar

% figure(4)
% clf
% hold all
% surf(r*1e9,t*1e9,u,'edgecolor','interp','facecolor','none')
% xlabel('Promieñ r [nm]','interpreter','latex')
% ylabel('Czas t [ns]','interpreter','latex')
% zlabel('St\c{e}\.zenie $y_{CO}$(r,t)[-]','interpreter','latex')
% grid on
% grid minor
% box on
% colormap(jet)
% colorbar
% view([-30 20])

for i=1:max(size(t))
    
    Rav(i)=3*trapz(r,k1.*u(i,:,1).*r.^2)/(dpm/2)^3;
    rav(i)=trapz(r,4.*gamma.*k1.*u(i,:,1).*r.^2);
    
end

figure(5)
clf
hold all
plot(r*1e9,u(end,:,3)*k1*gamma,'b-')
grid on
grid minor
box on
xlabel('Promieñ [nm]')
ylabel('Lokalne wartoœci szybkoœci reakcji odniesione do jednostki objêtoœci platyny')

% figure(6)
% clf
% hold all
% plot(t*1e9,rav,'c-','LineWidth',1.5)
% grid on
% grid minor
% box on
% xlabel('Czas t [ns]','interpreter','latex')
% ylabel('Szybko\''s\''c konwersji $\overline{r}_{CO}(t) \left[\frac{mol}{s}\right] $','interpreter','latex')
% 


function [c,f,s] = pdefun(r,t,u,dudr)
global k1 gamma Deff em
c = em/Deff;
f = dudr;
s = -gamma*k1/Deff*u;

function ic=IC(x)
ys0=0.001;

ic=ys0;

function [pl,ql,pr,qr] = BC(xl,ul,xr,ur,t)
ys0=0.001;
pl = 0;
ql = 1;
pr = ur - ys0;
qr = 0;

