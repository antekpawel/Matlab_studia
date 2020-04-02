clear 
clc
close all

Wczytny_obrz_z_pliku_tif = imread('1.tif');
%Uuwnie ottniej czeci mcierzy
Wczytny_obrz_z_pliku_tif(:,:,4) = [];
Konwerj_do_koloru_zrego=rgb2gray(Wczytny_obrz_z_pliku_tif);

% imshow(B)

Mp_pixeli_poidjaych_jnoc_powyzej_zdnej = Konwerj_do_koloru_zrego >= 150;

% imshow(Mp_pixeli_poidjaych_jnoc_powyzej_zdnej)

[y,x]=size(Mp_pixeli_poidjaych_jnoc_powyzej_zdnej);

for i=1:x
    for j=1:y
        if Mp_pixeli_poidjaych_jnoc_powyzej_zdnej(j,i)==1
            k(i)=j;
        end
    end
end
hold all
plot(k)

for i=1:x
    for j=y:(-1):1
        if Mp_pixeli_poidjaych_jnoc_powyzej_zdnej(j,i)==1
            g(i)=j;
        end
    end
end
plot(g)

% axis([0 1487 890 910])

