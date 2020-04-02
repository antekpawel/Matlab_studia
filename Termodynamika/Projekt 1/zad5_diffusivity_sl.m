function d=zad5_diffusivity_sl(v1,t,mi)

v1=v1.*1e6;
mi=mi.*1000;

d=2.98.*1e-11.*t.*mi.^(-1.026).*v1.^(-0.5473);