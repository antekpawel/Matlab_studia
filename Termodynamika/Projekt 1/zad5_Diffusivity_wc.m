function d=zad5_Diffusivity_wc(t,x,m,mi,v)

d=2.34.*1e-13.*t.*sqrt(x.*m)./mi./(v.*1e6).^0.6;