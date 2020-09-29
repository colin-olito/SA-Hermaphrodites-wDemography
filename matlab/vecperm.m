
function p = vecperm(m,n)
% vecperm
% function to calculate the vec permutation matrix of size m,n
% 4/9/03
%function p = vecperm(m,n)

p = zeros(m*n);
a = zeros(m,n);
for i = 1:m
    for j = 1:n
        e = a;
        e(i,j) = 1;
        p = p + kron(e,e');
    end
end


        
