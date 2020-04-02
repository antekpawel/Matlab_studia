function mini=optimummmm
dz=0.14*1.075;
di=dz:1e-4:0.4;
k=koszt(di);
df(di);

mini=fsolve(@(di) df(di),0.25);
Ti=Tizo(mini);
disp(mini)

function k=koszt(di)
di=0:0.001:0.5;
Z=1507;
dw=0.14;
dz=1.075.*dw;
ke=koszt_e(di);
ki=koszt_izo(di,Z,dz);
k=ke+ki;

function dff=df(di)

h=1e-6;
dff=(koszt(di-h)-koszt(di+h))./(2.*h);