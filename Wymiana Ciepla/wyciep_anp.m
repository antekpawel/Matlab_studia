function anp=wyciep_anp(Bi,minp)

pom=max(size(minp));

anp=zeros(pom,1);

for i=1:pom
    
    %anp(i)=(-1).^(i+1).*((2.*Bi.*sqrt(Bi.^2+minp(i).^2))./(minp(i).*(Bi.^2+Bi+minp(i).^2)));
    anp(i)=2.*sin(minp(i))./(minp(i)+sin(minp(i)).*cos(minp(i)));
    
end