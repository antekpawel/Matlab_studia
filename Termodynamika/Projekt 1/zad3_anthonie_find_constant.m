function [a,b,c,error,res]=zad3_anthonie_find_constant(p,t,prediction1,prediction2,prediction3)

%Relationship beetween vapour pressure and temperature lnPv=A+B/(T+C)
y = log(p);

%Anthonie function
F = @(x,t)x(1)+x(2)./(t+x(3));

%Predictions
x0 = [prediction1 prediction2 prediction3];

%Main
[x,error,res] = lsqcurvefit(F,x0,t,y);

%Constant
a=x(1);
b=x(2);
c=x(3);

