function d=zad5_diffusivity_hm(vk,t,mi)

v=vk*1e6;

eps=9.58./v-1.12;

mi=mi.*1000;

d=1.25.*1e-12.*(v.^(-0.19)-0.292).*t.^1.52.*mi.^eps;