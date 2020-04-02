function l=zad3_vapour_heat_Carruth(t,tc,omega)

R=8.3144;

%Reduced temperature
tr=t./tc;

%Equation
l=R.*tc.*(7.08.*(1-tr).^0.354+10.95.*omega.*(1-tr).^0.456);