function [x_bigA,x_a,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf,sm,delta)
                              
s1=theta(1,:);
s2=theta(2,:);
gam=theta(3,:);
lambda=zeros(3,1);
%selfing rate
C=0*ones(1,3);

    seldiff=[1 (1 - hf*sf) (1 - sf)];
    f=theta(4,:).*seldiff;
    seldiff_prime=[(1 - sm) (1 - hm*sm) 1];
    f_prime=theta(4,:).*seldiff_prime;

    for i=1:3
        U{i}=[s1(i)*(1-gam(i)) 0 ;
        s1(i)*gam(i) s2(i)];
        F{i}=(1-C(i))*[0 f(i);
        0 0];
        F_S{i}=C(i)*(1-delta)*[0 f(i);
        0 0];
        F_prime{i}=[0  f_prime(i)
        0 0];
    %Genotype specific population growth rate
        lambda(i)=max(eig(F{i}+U{i}+F_S{i}));
    end % for i
%Calculating the largest eigenvalue of the Jacobian for both boundaryies
Jac_bigA=[U{2}+.5.*F{2}+.5.*F_S{2}+.5.*(f_prime(2)/f_prime(1)).*F{1} F{3}+(f_prime(3)/f_prime(1)).*F{1}
            .25.*F_S{2} U{3}+F_S{3}];
Jac_a=[U{2}+.5.*F{2}+.5.*F_S{2}+.5*(f_prime(2)/f_prime(3)).*F{3} F{1}+(f_prime(1)/f_prime(3)).*F{3}
            .25.*F_S{2} U{1}+F_S{1}];    
%Dominant eigenvalue of the Jacobian at the AA boundary
x_bigA=max(eig(Jac_bigA))/lambda(1);
%Dominant eigenvalue of the Jacobian at the aa boundary
x_a=max(eig(Jac_a))/lambda(3);

lambda_bigA=lambda(1);
lambda_heterozyg=lambda(2);
lambda_a=lambda(3);

end