function alpha=zad1_alpha(C,lamp,d,U,vp,n)

alpha=C.*lamp./d.*(U.*d./vp).^n;