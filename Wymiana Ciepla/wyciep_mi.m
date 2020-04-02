function tab=wyciep_mi(Bi,how)

%0=ctg(mi)-mi/Bi

tab=zeros(1,how);

y=@(mi)tan(mi)-Bi/mi;

for i=1:how
    tab(i)=fzero(y,pi/2+pi.*i-pi);
end