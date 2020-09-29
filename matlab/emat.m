function [E] = emat(m,n,i,j)
%emat(m,n,i,j)
%matrix mXn, with 1 in i,j entry and zeros elsewhere

E=zeros(m,n);
E(i,j)=1;


end

