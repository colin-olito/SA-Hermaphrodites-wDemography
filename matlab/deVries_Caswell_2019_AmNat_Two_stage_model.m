%Program to implement Mendelian Matrix Population Models
%for two-stage model with parameters leading to evolutioary rescue
% Authors: Hal Caswell and Charlotte de Vries
% June 2018
clearvars; 
%Dimensions, 2 stages, 3 genotypes
om=2;
g=3;

%Identity and ones matrices
Ig=eye(g);
Iom=eye(om);
eg=ones(g,1);
eom=ones(om,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for EVOLUTIONARY SUICIDE (Fig 1A and 1C). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nzero=[0.1000         0   99.9000
    0.1000         0   99.9000];
%Baseline parameters
%theta=[s1, s2, g, f]

theta=[0.65 0.7 0.05 3.2]'*ones(1,3);
%Selection differentials

% %s1, juvenile survival
ds1=0.032;
 seldiff=[-ds1+0.01 ds1+0.008 -ds1+0.005];
 theta(1,:)=theta(1,:)+seldiff;

%s2 differential, adult survival
% s2
    ds2=0.0001;

seldiff=[0 0 ds2];
theta(2,:)=theta(2,:)+seldiff;

%this is the f differential for the pleiotropy example
% f
    df=0.81;
    seldiff=[df/2 -df df];
theta(4,:)=theta(4,:)+seldiff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% END of Parameters for EVOLUTIONARY SUICIDE (Fig 1A and 1C). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for EVOLUTIONARY RESCUE (Fig 1B and 1D). 
% Uncommment the following section of code if you wish to use the
% parameters from Figure 1B and 1D

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p=0.03;
% %Initial condition
% Nzero=[p 0 1-p
%        p 0 1-p];
% %Baseline parameters
% %theta=[s1, s2, g, f]
% theta=[0.7 0.7 0.05 3.2]'*ones(1,3);
% %Selection Differentials, juvenile survival
% % %s1 
%  seldiff=[0.032 -0.04 0.027];
%  theta(1,:)=theta(1,:)+seldiff;
% 
% % s2, adult survival
%     ds4=0.0001;
% seldiff2=[0 0 -ds4];
% theta(2,:)=theta(2,:)+seldiff2;
% %Fertility differential
% % f
%     df=0.81;
%     seldiff=[-df/2 df -df];
% theta(4,:)=theta(4,:)+seldiff;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% END of Parameters for EVOLUTIONARY RESCUE (Fig 1B and 1D). 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Simulation time
tlimit=4000;



%Parameter assignment 
s1=theta(1,:);
s2=theta(2,:);
gam=theta(3,:);
f=theta(4,:);
lambda=zeros(3,1);

for i=1:3
    U{i}=[s1(i)*(1-gam(i)) 0 ;
        s1(i)*gam(i) s2(i)];
    F{i}=[0 f(i);
        0 0];
    %Genotype specific population growth rate
    lambda(i)=max(eig(F{i}+U{i}))
end % for i


%Calculating the largest eigenvalue of the Jacobian for both boundaryies

inv_bigA=max(eig((U{2}+0.5*F{1}+0.5*F{2}))); 
inv_a=max(eig(U{2}+0.5*F{3}+0.5*F{2}));

%Dominant eigenvalue of the Jacobian at the AA boundary
x_bigA=inv_bigA/lambda(1)
%Dominant eigenvalue of the Jacobian at the aa boundary
x_a=inv_a/lambda(3)

d=eye(g);
c=[0; 1]; %Indicator vector, with a 1 in positions corresponding to breeding stages

bbD=kron(Iom,d); %Equation 27 main text
bbU=zeros(g*om); 
bbF=zeros(g*om);
for i=1:g
    bbU=bbU+kron(emat(g,g,i,i),U{i}); %Equation 25 main text
    bbF=bbF+kron(emat(g,g,i,i),F{i}); %Equation 26 main text
    
end % for i

%Matrices
Ig=eye(g);
K=vecperm(om,g);

X=kron(Ig,c'); %Equation 18

W=[1 .5 0
    0 .5 1];

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
nzero=vec(Nzero);
n=nzero;
    
    for i=1:tlimit
        
        nout(:,i)=n;
                
        %breeding population values
        Nb=sum(X*n);
        pb=X*n/Nb;
        qb=W*pb;
        rb=M*qb;
        
        Utilde=K'*bbD*K*bbU;
        
        h=zeros(g);
        for i=1:g
            pi=Ig(:,i);
            qi=W*pi;
            ri=M*qi;
            
            piprime=Z*kron(ri,rb);
            
            H(:,i)=piprime;
                       
        end
        bbH=kron(Iom,H);
       
        Ftilde=K'*bbH*K*bbF;
        nprime=(Utilde + Ftilde)*n;
        pprime=X*nprime/norm(X*nprime,1);
        pout=[pprime];
        
        n=nprime;
        p=pprime;
        
    end
    
    
    
    
    %change output name
    n=nout;
    clear nout;
    
nsum=sum(n);
for i=1:4000;
pfull(:,i)=n(:,i)./nsum(i);
end



%Figures

%set(0,'defaultaxeslinewidth','factory');
set(0,'defaultaxesfontsize',18);
%set(0,'defaultaxeslinestyle','factory')
set(0,'defaultaxescolororder','factory');

set(0,'defaultaxescolororder',[1 0 0; 0 1 0; 0 0 1;0 0 0])
% set(0,'DefaultAxesLinestyleOrder','-|--|:|-o|-*|-v|-^');
set(0,'DefaultLineLineWidth',2)
set(0,'DefaultLineMarkerSize',2)

subplot(2,1,1);
semilogy(1:tlimit,kron(eg',Iom)*n); %Population abundance in each stage, summed over genotypes
xlabel('Time')
ylabel('Abundance')
legend('Juvenile','Adult','location','best')
%title('Evolutionary rescue')

b=kron(Ig,eom')*n;
pstar=b*diag(1./sum(n));

ntemp=W*kron(Ig,eom')*n;
qstar=ntemp*diag(1./sum(ntemp));

subplot(2,1,2);
plot(1:tlimit,pstar) %Genotype frequencies, aggregated over all stages
xlabel('Time')
ylabel('Frequency')
legend('AA','Aa','aa','location','best')
title('$\lambda_{Aa}>\lambda_{AA}$ and $\lambda_{Aa}>\lambda_{aa}$','interpreter','latex');
