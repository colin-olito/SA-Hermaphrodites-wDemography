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

%Params from Colins R code: fwdDemModelSim  <-  function(	om = 2, g = 3, theta = c(0.6,0.6,0.05,5.9), theta_prime = c(0.6,0.6,0.05,5.9), 
%								hf = 1/2, hm = 1/2, sf = 0.1, sm = 0.105, C = 0, delta = 0, tlimit = 10^4, Ainvade = FALSE, intInit = FALSE, ...) {

Nzero=[99.9000 0 0.1000
    99.9000 0 0.1000];
%Baseline parameters
%theta=[s1, s2, g, f]

theta=[0.6 0.6 0.05 6.1]'*ones(1,3);
sm=0:0.0001:0.15;
lambda_save_bigA=nan*ones(length(sm),3);
lambda_save_a=nan*ones(length(sm),3);

%The precision below which the algorithm stops searching for the boundary sm
%value
deltathres=0.00001;

%Selection differentials
%this is the f differential for the pleiotropy example
% f
    hf=.25;
    hm=.25;
    
Thres_bigA=nan*ones(length(sm),1);
Thres_a=nan*ones(length(sm),1);

delta=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find threshold zeta_AA=1 for each sm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:length(sm)

    if isnan(Thres_bigA(j)) 
        sf_temp=[0 0.15]; delta_titrate=1;
                    % first just a  check
                    [x_bigA0,x_a0,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,0,sm(j),delta);
                    [x_bigA1,x_a1,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,0.15,sm(j),delta);% who wins when d=0 (one would think fast does, unless slowfecfactor is unreasonably high)
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
                    elseif x_bigA0>=1  && x_bigA1<=1 % there will be a threshold, let's look for it
                        while delta_titrate>deltathres
                            delta_titrate=delta_titrate/2;
                                [x_bigA,x_a,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,mean(sf_temp),sm(j),delta);                            
                            if x_bigA>1 % big A unstable at halfway point, so true threshold must be somewhere between here and the current upper
                                sf_temp=[mean(sf_temp) sf_temp(2)];
                            else % slow won at halfway point, so true threshold must be somewhere between here and the current lower
                                sf_temp=[sf_temp(1) mean(sf_temp)];
                            end
                        end
                     mean(sf_temp);
                    Thres_bigA(j)=mean(sf_temp);
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
for j=1:length(sm)

    if isnan(Thres_a(j)) 
        sf_temp=[0 0.15]; delta_titrate=1;
                    % first just a  check
                    [x_bigA0,x_a0,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,0,sm(j),delta);
                    [x_bigA1,x_a1,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,0.15,sm(j),delta);% who wins when d=0 (one would think fast does, unless slowfecfactor is unreasonably high)
                    if x_a0<1 && x_a1<1 % a wins across the entire range of sm values; let's mark this with an above-range value of 1.1
                         Thres_bigA(j)=nan;
                    elseif x_a0>1 && x_a1>1 % a unstable across the entire range of sm values; let's mark this with a below-range value of -0.1
                        Thres_bigA(j)=nan;
                    elseif x_a1>=1  && x_a0<=1 % there will be a threshold, let's look for it
                        while delta_titrate>deltathres
                            delta_titrate=delta_titrate/2;
                                [x_bigA,x_a,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,mean(sf_temp),sm(j),delta);                            
                            if x_a>1 % a unstable at halfway point, so true threshold must be somewhere between here and the current lower
                                sf_temp=[sf_temp(1) mean(sf_temp)];
                            else % a stable at halfway point, so true threshold must be somewhere between here and the current lower
                                sf_temp=[mean(sf_temp) sf_temp(2)];
                            end
                        end
                    Thres_a(j)=mean(sf_temp);
                    lambda_save_a(j,1)=lambda_bigA;
                    lambda_save_a(j,2)=lambda_heterozyg;
                    lambda_save_a(j,3)=lambda_a;
                    else
                        disp('something is odd'); pause
                    end


    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find whether eqm is exponentially increasing or decreasing for each sm,
% sf combination

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nzero=[1 1 1 1 1 1]';
sf=0:0.0001:0.15;
lambda_eq=nan*ones(length(sm),length(sf));
for i=1:length(sm)
    for j=1:length(sf)
        if sf(j)<Thres_a(i)
         [x_bigA,x_a,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),sm(i),delta);
        lambda_eq(j,i)=lambda_a;
        elseif sf(j)>Thres_bigA(i)
        [x_bigA,x_a,lambda_bigA,lambda_heterozyg,lambda_a]=eigenvalues_Jacobian(theta,hf,hm,sf(j),sm(i),delta);
        lambda_eq(j,i)=lambda_bigA;
        else
             [lambda_sim]=simulate_dyn(theta,hf,hm,sf(j),sm(i),delta,nzero);
             lambda_eq(j,i)=lambda_sim;
        end  
    end
end

