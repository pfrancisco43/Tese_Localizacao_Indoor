function [estEM,estES,erroEM,erroES]=main(BS,SP,UE,sigmaT,sigmaA)
    %replicação do artigo A Simple Method for 5G Positioning and 
    %Synchronization without Line-of-Sight. 
    %Autor: WYMEERSCH, H. 
    
    sigma.TOA = sigmaT; 
    sigma.AOA_az = sigmaA;    %[radian]  standard deviation of DOA in azimuth
    sigma.AOA_el = sigmaA;    %[radian]  standard deviation of DOA in elevation
    sigma.AOD_az = sigmaA;    %[radian]  standard deviation of DOD in azimuth
    sigma.AOD_el = sigmaA;    %[radian]  standard deviation of DOD in elevation
    alfa=pi/3+0.06;             %direction
    B=20+0.8;                   %Bias
    covariancia=diag([sigma.TOA^2,sigma.AOD_az^2,sigma.AOD_el^2,sigma.AOA_az^2,sigma.AOA_el^2]);
    [~,medidas]=geraMedidas(BS,UE,SP,covariancia,B,alfa);
    
       
    pos=BS;
    BS=[];
    BS.pos=pos;
    pos=UE;
    UE=[];
    UE.state=pos;
    UE.state(4)=alfa;
    UE.state(5)=B;
    
    SPO=SP;
    pos=SP;
    SP=[];
    N_SP=size(pos,2);
    for k=1:N_SP
        SP(k).pos=pos(:,k);
    end
    
    % 0. Simulation parameters
    % ------------------------
    
    bias_prior=200.00;   % [m] uncertainty of the bias (-bias_prior/2,+bias_prior/2)
    alpha_prior=pi;      % [rad] uncertainty of the orientation (0,alpha_prior);
    Ns=5;               % number of samples in the objective function
    Na=20;               % grid points for user orientation in [0, pi]
    
    Nb=5;
    
    % 2. Determine an estimate for bias and orientation
    % --------------------------------------------------
    f = @(x)AlgoritmoHenk.computeErrorMetric(x,Ns,covariancia,medidas,BS);

    % 3. For the chosen bias and orientation, determine UE position
    % --------------------------------------------------------------
    
    xbest=[20 pi/3];
    [dum,UE.mean,UE.cov] = f(xbest);
        
    % 4. Determine the SP locations
    % ------------------------------
    SP=AlgoritmoHenk.mapEnvironment(BS,SP,UE.mean,UE.cov,medidas);
        
    erroEM=norm(UE.mean(1:3)-UE.state(1:3));
    estEM=UE.mean(1:3);
            
    for i=1:N_SP
        si(i)=norm(SP(i).pos - SP(i).mean);
        estES(:,i)=SP(i).mean;
    end
    erroES=mean(si);