%Program to implement Mendelian Matrix Population Models
%for two-stage model with selfing
%Version where how-you-were-produced does not affect your demographic rates
%Authors: C. de Vries and C. Olito
% July 2020
clearvars ;
%Dimensions, 2 stages, 3 genotypes
om=2;
g=3;

%Identity and ones matrices
Ig=eye(g);
Iom=eye(om);
eg=ones(g,1);
eom=ones(om,1);

%Params from Colins R code: fwdDemModelSim  <-  function(	om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = c(0.6,0.6,0.05,5.9), 
%								hf = 1/2, hm = 1/2, sf = 0.1, sm = 0.105, C = 0, delta = 0, tlimit = 10^4, Ainvade = FALSE, intInit = FALSE, ...) {

Nzero=[99.9000 0 0.1000
    99.9000 0 0.1000];
%Baseline parameters
%theta=[s1, s2, g, f]

theta=[0.6 0.6 0.05 5.8]'*ones(1,3);
sf=0:0.001:0.15;
lambda_save_bigA=nan*ones(length(sf),3);
lambda_save_a=nan*ones(length(sf),3);

%The precision below which the algorithm stops searching for the boundary sm
%value
deltathres=0.0001;

%Selection differentials
%this is the f differential for the pleiotropy example
% f
    hf=.25;
    hm=.25;
    
Thres_bigA=nan*ones(length(sf),1);
Thres_a=nan*ones(length(sf),1);

delta=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find threshold zeta_AA=1 for each sm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(sf)

    if isnan(Thres_bigA(j)) 
        sm_temp=[0 0.15]; delta_titrate=1;
                    % first just a  check
                    [x_bigA0,x_a0,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),0,delta);
                    [x_bigA1,x_a1,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),0.15,delta);
                    if x_bigA0<1 && x_bigA1<1 % big A wins across the entire range of sm values; let's mark this with an above-range value of 1.1
                         Thres_bigA(j)=nan;
                    lambda_save_bigA(j,1)=lambda_bigA;
                    lambda_save_bigA(j,2)=lambda_heterozyg;
                    lambda_save_bigA(j,3)=lambda_a;
                    elseif x_bigA0>1 && x_bigA1>1 % big A unstable across the entire range of sm values; let's mark this with a below-range value of -0.1
                        Thres_bigA(j)=nan;
                    lambda_save_bigA(j,1)=lambda_bigA;
                    lambda_save_bigA(j,2)=lambda_heterozyg;
                    lambda_save_bigA(j,3)=lambda_a;
                    elseif x_bigA0<=1  && x_bigA1>=1 % there will be a threshold, let's look for it
                        while delta_titrate>deltathres
                            delta_titrate=delta_titrate/2;
                                [x_bigA,x_a,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),mean(sm_temp),delta);                            
                            if x_bigA>1 % big A unstable at halfway point, so true threshold must be somewhere between here and the current lower
                                sm_temp=[sm_temp(1) mean(sm_temp)];
                            else % slow won at halfway point, so true threshold must be somewhere between here and the current lower
                                sm_temp=[mean(sm_temp) sm_temp(2)];
                            end
                        end
                    Thres_bigA(j)=mean(sm_temp);
                    lambda_save_bigA(j,1)=lambda_bigA;
                    lambda_save_bigA(j,2)=lambda_heterozyg;
                    lambda_save_bigA(j,3)=lambda_a;
                    else
                        disp('something is odd'); pause
                    end


    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find threshold zeta_aa=1 for each sm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=1:length(sf)

    if isnan(Thres_a(j)) 
        sm_temp=[0 0.15]; delta_titrate=1;
                    % first just a  check
                    [x_bigA0,x_a0,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),0,delta);
                    [x_bigA1,x_a1,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),0.15,delta);
                    if x_a0<1 && x_a1<1 % a wins across the entire range of sm values; let's mark this with an above-range value of 1.1
                         Thres_a(j)=nan;
                    elseif x_a0>1 && x_a1>1 % a unstable across the entire range of sm values; let's mark this with a below-range value of -0.1
                        Thres_a(j)=nan;
                    elseif  x_a0>=1  && x_a1<=1  % there will be a threshold, let's look for it
                        while delta_titrate>deltathres
                            delta_titrate=delta_titrate/2;
                                [x_bigA,x_a,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),mean(sm_temp),delta);                            
                            if x_a>1 % a unstable at halfway point, so true threshold must be somewhere between here and the current upper
                                sm_temp=[mean(sm_temp) sm_temp(2)];
                            else % a stable at halfway point, so true threshold must be somewhere between here and the current lower
                                sm_temp=[sm_temp(1) mean(sm_temp)];
                            end
                        end
                    Thres_a(j)=mean(sm_temp);
                    lambda_save_a(j,1)=lambda_bigA;
                    lambda_save_a(j,2)=lambda_heterozyg;
                    lambda_save_a(j,3)=lambda_a;
                    else
                        disp('something is odd'); pause
                    end


    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find extinction boundaries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find extinction boundary inside coexistence region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thres_extinct_coex=nan*ones(length(sm),1);
% sf_upper=nan*ones(length(sm),1);
% options = optimoptions('fsolve');
% nzero=[1 1 1 1 1 1]';
% for j=1:length(sm)
% 
%     if isnan(Thres_extinct_coex(j)) 
%         if isnan(Thres_bigA(j)) 
%             sf_upper(j)=0.15;
%         else
%             sf_upper(j)=Thres_bigA(j); 
%         end
%        Thres_extinct_coex(j) = patternsearch(@(x)simulate_dyn_solve(x,sm(j)),sm(j),[],[],[],[],Thres_a(j),sf_upper(j),options);
%       
%     end
% 
% end

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find extinction boundary inside coexistence region

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Thres_extinct_coex_flip=nan*ones(length(sf),1);
sm_upper=nan*ones(length(sf),1);
sf=0:0.001:0.15;

nzero=[1 1 1 1 1 1]';
for j=1:length(sf)

    if isnan(Thres_extinct_coex_flip(j)) 
        if isnan(Thres_a(j))
            sm_upper(j)=0.15;
        else
            sm_upper(j)=Thres_a(j); 
        end
         sm_temp=[Thres_bigA(j) sm_upper(j)];
          delta_titrate=1;
                    % first just a  check
                   [lambda_sim0]=simulate_dyn(theta,hf,hm,sf(j),sm_temp(1),delta,nzero);
                   [lambda_sim1]=simulate_dyn(theta,hf,hm,sf(j),sm_temp(2),delta,nzero);
                    if lambda_sim0>1 && lambda_sim1>1 % No extinction, positive growth rate for whole range
                         Thres_extinct_coex_flip(j)=nan;
                    elseif  lambda_sim0<1 && lambda_sim1<1  % extinction across the entire range of sm values; let's mark this with a below-range value of -0.1
                        Thres_extinct_coex_flip(j)=0;
                    else  % there is a threshold, let's look for it
                        while delta_titrate>deltathres
                            delta_titrate=delta_titrate/2;
                                    lambda_sim=simulate_dyn(theta,hf,hm,sf(j),mean(sm_temp),delta,nzero);
                            if lambda_sim>1 % viable at halfway point, so true threshold must be somewhere between here and the current upper value
                                sm_temp=[mean(sm_temp) sm_temp(2)];
                            else % a stable at halfway point, so true threshold must be somewhere between here and the current lower
                                sm_temp=[sm_temp(1) mean(sm_temp)];
                            end
                        end
                    Thres_extinct_coex_flip(j)=mean(sm_temp);
%                     else
%                         disp('something is odd'); pause
                    end
  

    end

end

% for j=1:length(sm)
% 
%     if isnan(Thres_extinct_coex(j)) 
%         if isnan(Thres_bigA(j))
%             sf_upper(j)=0.15;
%         else
%             sf_upper(j)=Thres_bigA(j); 
%         end
%          sf_temp=[Thres_a(j) sf_upper(j)];
%           delta_titrate=1;
%                     % first just a  check
%                    [lambda_sim0]=simulate_dyn(theta,hf,hm,sf_temp(1),sm(j),delta,nzero);
%                    [lambda_sim1]=simulate_dyn(theta,hf,hm,sf_temp(2),sm(j),delta,nzero);
%                     if lambda_sim0>=1 && lambda_sim1>=1 % No extinction, positive growth rate for whole range
%                          Thres_extinct_coex(j)=0.15;
%                     elseif  lambda_sim0<1 && lambda_sim1<1  % extinction across the entire range of sm values; let's mark this with a below-range value of -0.1
%                         Thres_extinct_coex(j)=0;
%                     else  % there is a threshold, let's look for it
%                         while delta_titrate>deltathres
%                             delta_titrate=delta_titrate/2;
%                                     lambda_sim=simulate_dyn(theta,hf,hm,mean(sf_temp),sm(j),delta,nzero);
%                             if lambda_sim>1 % viable at halfway point, so true threshold must be somewhere between here and the current upper value
%                                 sf_temp=[sf_temp(1) mean(sf_temp)];
%                             else % a stable at halfway point, so true threshold must be somewhere between here and the current lower
%                                 sf_temp=[mean(sf_temp) sf_temp(2)];
%                             end
%                         end
%                     Thres_extinct_coex(j)=mean(sf_temp);
% %                     else
% %                         disp('something is odd'); pause
%                     end
%   
% 
%     end
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find extinction boundary inside aa region 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Thres_extinct_aa=nan*ones(length(sf),1);
sm_upper=nan*ones(length(sf),1);
nzero=[1 1 1 1 1 1]';
for j=1:length(sf)

    if isnan(Thres_extinct_aa(j)) 

         sm_temp=[Thres_a(j) 0.15];
          delta_titrate=1;
                    % first just a  check
                   [lambda_sim0]=simulate_dyn(theta,hf,hm,sf(j),sm_temp(1),delta,nzero);
                   [lambda_sim1]=simulate_dyn(theta,hf,hm,sf(j),sm_temp(2),delta,nzero);
                    if lambda_sim0>=1 && lambda_sim1>=1 % No extinction, positive growth rate for whole range
                         Thres_extinct_aa(j)=.15;
                    elseif  lambda_sim0<1 && lambda_sim1<1  % extinction across the entire range of sm values; let's mark this with a below-range value of -0.1
                        Thres_extinct_aa(j)=Thres_a(j);
                    elseif lambda_sim0>=1 && lambda_sim1<1  % there is a threshold, let's look for it
                        while delta_titrate>deltathres
                            delta_titrate=delta_titrate/2;
                                    lambda_sim=simulate_dyn(theta,hf,hm,sf(j),mean(sm_temp),delta,nzero);
                            if lambda_sim>1 % viable at halfway point, so true threshold must be somewhere between here and the current upper value
                                sm_temp=[mean(sm_temp) sm_temp(2)];
                            else % a stable at halfway point, so true threshold must be somewhere between here and the current lower
                                sm_temp=[sm_temp(1) mean(sm_temp)];
                            end
                        end
                    Thres_extinct_aa(j)=mean(sm_temp);
                    else
                        disp('something is odd'); pause
                    end
  

    end

end


 save titrate_attempt1_flip

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find whether eqm is exponentially increasing or decreasing for each sm,
% sf combination

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% nzero=[1 1 1 1 1 1]';
% sf=0:0.001:0.15;
% lambda_eq=nan*ones(length(sm),length(sf));
% for i=1:length(sm)
%     for j=1:length(sf)
%         if sf(j)<Thres_a(i)
%          [x_bigA,x_a,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),sm(i),delta);
%         lambda_eq(j,i)=lambda_a;
%         elseif sf(j)>Thres_bigA(i)
%         [x_bigA,x_a,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),sm(i),delta);
%         lambda_eq(j,i)=lambda_bigA;
%         else
%              [lambda_sim]=simulate_dyn(theta,hf,hm,sf(j),sm(i),delta,nzero);
%              lambda_eq(j,i)=lambda_sim;
%         end  
%     end
% end
% 
%  save titrate_attempt1
