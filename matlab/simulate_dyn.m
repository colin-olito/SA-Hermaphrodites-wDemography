function [lambda_sim]=simulate_dyn(theta, hf,hm,sf,sm,delta,nzero)
 %[lambda_sim,pout,nout]=simulate_dyn(theta,hf,hm,sf,sm,delta,nzero)
                      
s1=theta(1,:);
s2=theta(2,:);
gam=theta(3,:);
lambda=zeros(3,1);
%selfing rate
C=0*ones(1,3);
g=3;
om=2;
Ig=eye(g);
Iom=eye(om);
eg=ones(g,1);
eom=ones(om,1);

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


d=eye(g);
c=[0; 1]; %Indicator vector, with a 1 in positions corresponding to breeding stages

bbD=kron(Iom,d); %Equation 27 main text
bbU=zeros(g*om); 
bbF=zeros(g*om);
bbF_S=zeros(g*om);
bbF_prime=zeros(g*om);
for i=1:g
    bbU=bbU+kron(emat(g,g,i,i),U{i}); %Equation 25 main text
    bbF=bbF+kron(emat(g,g,i,i),F{i}); %Equation 26 main text
    bbF_prime=bbF_prime+kron(emat(g,g,i,i),F_prime{i});
    bbF_S=bbF_S+kron(emat(g,g,i,i),F_S{i});
end % for i

%Matrices
Ig=eye(g);
K=vecperm(om,g);

X=kron(Ig,c'); %Equation 18

W=[1 .5 0
    0 .5 1];

W_prime=[ones(1,om) .5*ones(1,om) 0*ones(1,om)
    0*ones(1,om) .5*ones(1,om) ones(1,om)];

Z=[1 0 0 0
    0 1 1 0
    0 0 0 1];

%No mutations
u=0;
v=0;
M=[1-u v;
    u 1-v];

bbM=kron(Ig,M);

Q=Z*kron(W,W)*kron(Ig,X); %Equation 21
n=nzero;
tlimit=2000;
    
nout=nan*ones(tlimit,length(nzero)); 
pout=nan*ones(tlimit,length(nzero));
it=0;
thres=0.000001;
pdiff=1;
while it<tlimit && pdiff>thres
        it=it+1;
        nout(it,:)=n;
                
        %breeding population values
%         Nb=sum(X*n);
%         pb=X*n/Nb;
%         qb=W*pb;
%         rb=M*qb;
        
        Utilde=K'*bbD*K*bbU;
        ngam=ones(1,2)*W_prime*bbF_prime*n;
        q_prime=W_prime*bbF_prime*n/ngam;
        rb=M*q_prime;
        
        h=zeros(g);
        for i=1:g
            pi=Ig(:,i);
            qi=W*pi;
            ri=M*qi;
            
            piprime=Z*kron(ri,rb);
            
            H(:,i)=piprime;
                       
        end
        
        H_S=[1 .25 0 
            0 .5 0
            0 .25 1];
        
        bbH=kron(Iom,H);
        bbH_S=kron(Iom,H_S);
        Ftilde=K'*bbH*K*bbF;
        Ftilde_selfing=K'*bbH_S*K*bbF_S;
        nprime=(Utilde + Ftilde+Ftilde_selfing)*n;
        n=nprime;
        pout(it,:)=n./sum(n);
       
       % nout(it,:)=n;
        if it>3
        pdiff=sqrt(sum((pout(it,:) - pout(it-1,:)) .^ 2));
        lambda_sim=sum(nout(it,:))./sum(nout(it-1,:));
        end
end

end