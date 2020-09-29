%Program to implement Mendelian Matrix Population Models
%for two-stage model with selfing
%Version where how-you-were-produced does not affect your demographic rates
%Authors: C. de Vries and C. Olito
% July 2020
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
%Params from Colins R code: fwdDemModelSim  <-  function(	om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = c(0.6,0.6,0.05,5.9), 
%								hf = 1/2, hm = 1/2, sf = 0.1, sm = 0.105, C = 0, delta = 0, tlimit = 10^4, Ainvade = FALSE, intInit = FALSE, ...) {

Nzero=[99.9000 0 0.1000
    99.9000 0 0.1000];
%Baseline parameters
%theta=[s1, s2, g, f]

theta=[0.6 0.6 0.05 6.1]'*ones(1,3);
%Selection differentials

% %s1, juvenile survival
% ds1=0.032;
%  seldiff=[-ds1+0.01 ds1+0.008 -ds1+0.005];
%  theta(1,:)=theta(1,:)+seldiff;

%s2 differential, adult survival
% s2
%     ds2=0.0001;
% 
% seldiff=[0 0 ds2];
% theta(2,:)=theta(2,:)+seldiff;

%this is the f differential for the pleiotropy example
% f
    hf=.5;
    hm=.5;
    sf=0.1;
    sm=0.105;
    seldiff=[1 (1 - hf*sf) (1 - sf)];
theta(4,:)=theta(4,:).*seldiff;

seldiff_prime=[(1 - sm) (1 - hm*sm) 1];
f_prime=theta(4,:).*seldiff_prime;
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
tlimit=7000;



%Parameter assignment 
s1=theta(1,:);
s2=theta(2,:);
gam=theta(3,:);
f=theta(4,:);
lambda=zeros(3,1);
%selfing rate
C=0.25*ones(1,3);

%selfing fitness reduction
maxi=10;
x_bigA=zeros(maxi,1);
x_a=zeros(maxi,1);
lambda_save=zeros(maxi,3);
delta_save=zeros(maxi,1);

for ii=1:maxi
delta=ii*(.9/maxi);
delta_save(ii,1)=delta;

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
    lambda_save(ii,i)=max(eig(F{i}+U{i}+F_S{i}));
end % for i


%Calculating the largest eigenvalue of the Jacobian for both boundaryies

% inv_bigA=max(eig((U{2}+0.5*F{1}+0.5*F{2}))); 
% inv_a=max(eig(U{2}+0.5*F{3}+0.5*F{2}));

Jac_bigA=[U{2}+.5.*F{2}+.5.*F_S{2}+.5.*(f_prime(2)/f_prime(1)).*F{1} F{3}+(f_prime(3)/f_prime(1)).*F{1}
            .25.*F_S{2} U{3}+F_S{3}];
Jac_a=[U{2}+.5.*F{2}+.5.*F_S{2}+.5*(f_prime(2)/f_prime(3)).*F{3} F{1}+(f_prime(1)/f_prime(3)).*F{3}
            .25.*F_S{2} U{1}+F_S{1}];

        
%Dominant eigenvalue of the Jacobian at the AA boundary
x_bigA(ii,1)=max(eig(Jac_bigA))/lambda_save(ii,1);
%Dominant eigenvalue of the Jacobian at the aa boundary
x_a(ii,1)=max(eig(Jac_a))/lambda_save(ii,3);


end
% 
% d=eye(g);
% c=[0; 1]; %Indicator vector, with a 1 in positions corresponding to breeding stages
% 
% bbD=kron(Iom,d); %Equation 27 main text
% bbU=zeros(g*om); 
% bbF=zeros(g*om);
% bbF_S=zeros(g*om);
% bbF_prime=zeros(g*om);
% for i=1:g
%     bbU=bbU+kron(emat(g,g,i,i),U{i}); %Equation 25 main text
%     bbF=bbF+kron(emat(g,g,i,i),F{i}); %Equation 26 main text
%     bbF_prime=bbF_prime+kron(emat(g,g,i,i),F_prime{i});
%     bbF_S=bbF_S+kron(emat(g,g,i,i),F_S{i});
% end % for i
% 
% %Matrices
% Ig=eye(g);
% K=vecperm(om,g);
% 
% X=kron(Ig,c'); %Equation 18
% 
% W=[1 .5 0
%     0 .5 1];
% 
% W_prime=[ones(1,om) .5*ones(1,om) 0*ones(1,om)
%     0*ones(1,om) .5*ones(1,om) ones(1,om)];
% 
% Z=[1 0 0 0
%     0 1 1 0
%     0 0 0 1];
% 
% %No mutations
% u=0;
% v=0;
% M=[1-u v;
%     u 1-v];
% 
% bbM=kron(Ig,M);
% 
% Q=Z*kron(W,W)*kron(Ig,X); %Equation 21
% nzero=vec(Nzero);
% n=nzero;
%     
%     for i=1:tlimit
%         
%         nout(:,i)=n;
%                 
%         %breeding population values
% %         Nb=sum(X*n); pb=X*n/Nb; qb=W*pb; rb=M*qb;
%         
%         Utilde=K'*bbD*K*bbU;
%         ngam=ones(1,2)*W_prime*bbF_prime*n;
%         q_prime=W_prime*bbF_prime*n/ngam;
%         rb=M*q_prime;
%         
%         h=zeros(g);
%         for i=1:g
%             pi=Ig(:,i);
%             qi=W*pi;
%             ri=M*qi;
%             
%             piprime=Z*kron(ri,rb);
%             
%             H(:,i)=piprime;
%                        
%         end
%         
%         H_S=[1 .25 0 
%             0 .5 0
%             0 .25 1];
%         
%         bbH=kron(Iom,H);
%         bbH_S=kron(Iom,H_S);
%         Ftilde=K'*bbH*K*bbF;
%         Ftilde_selfing=K'*bbH_S*K*bbF_S;
%         nprime=(Utilde + Ftilde+Ftilde_selfing)*n;
%         pprime=X*nprime/norm(X*nprime,1);
%         pout=[pprime];
%         
%         n=nprime;
%         p=pprime;
%         
%     end
%     
%     
%     
%     
%     %change output name
%     n=nout;
%     clear nout;
%     
% nsum=sum(n);
% for i=1:4000;
% pfull(:,i)=n(:,i)./nsum(i);
% end
% 
% 
% 
% %Figures
% 
% %set(0,'defaultaxeslinewidth','factory');
% set(0,'defaultaxesfontsize',18);
% %set(0,'defaultaxeslinestyle','factory')
% set(0,'defaultaxescolororder','factory');
% 
% set(0,'defaultaxescolororder',[1 0 0; 0 1 0; 0 0 1;0 0 0])
% % set(0,'DefaultAxesLinestyleOrder','-|--|:|-o|-*|-v|-^');
% set(0,'DefaultLineLineWidth',2)
% set(0,'DefaultLineMarkerSize',2)
% 
% subplot(2,1,1);
% semilogy(1:tlimit,kron(eg',Iom)*n); %Population abundance in each stage, summed over genotypes
% xlabel('Time')
% ylabel('Abundance')
% legend('Juvenile','Adult','location','best')
% %title('Evolutionary rescue')
% 
% b=kron(Ig,eom')*n;
% pstar=b*diag(1./sum(n));
% 
% ntemp=W*kron(Ig,eom')*n;
% qstar=ntemp*diag(1./sum(ntemp));
% 
% subplot(2,1,2);
% plot(1:tlimit,pstar) %Genotype frequencies, aggregated over all stages
% xlabel('Time')
% ylabel('Frequency')
% legend('AA','Aa','aa','location','best')
% title('$\lambda_{Aa}>\lambda_{AA}$ and $\lambda_{Aa}>\lambda_{aa}$','interpreter','latex');
