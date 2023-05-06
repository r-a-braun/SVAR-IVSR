clear; clc; addpath('functions')

%% Computations for Panel A and B
load('results\output_BH19_R1.mat') 
n = 4; ntil = 5;
nrep = output_R1.input.nrep; 
for i = 1:nrep
    B_HR19 = reshape(output_R1.vecBs(:,i),ntil,ntil); 
    Atil = inv(B_HR19);
    A = Atil(1:4,1:4);
    A = A./[A(1,1),A(2,2),A(3,1),A(4,4)]';
    alpha_qp_R1(i) = - A(1,3);    
end
load('results/output_BH19_R2.mat')
n = 4; ntil = 5;
nrep = output_BH19_R2.input.nrep;
probs_HR20 = zeros(2,nrep);
probs_BH19 = zeros(2,nrep);
c_alpha_qp = .1; sigma_alpha_qp = 0.1; nu_alpha_qp = 3;
for i = 1:nrep
    B_HR19 = reshape(output_BH19_R2.vecBs(:,i),ntil,ntil);
    B_HR19 = B_HR19(:,[n,1:n-1,n+1]);
    Atil = inv(B_HR19);
    A = Atil(1:4,1:4);
    A = A./[A(1,1),A(2,2),A(3,1),A(4,4)]';
    alpha_qp_R2(i) = - A(1,3);     
    probs_HR20(1,i) = alpha_qp_R2(i)<0.04 ;
    probs_BH19(1,i) = student_pos_prior(alpha_qp_R2(i) ,c_alpha_qp,sigma_alpha_qp,nu_alpha_qp);
    %PRIOR
    B_HR19 = reshape(output_BH19_R2.vecBs_pr(:,i),ntil,ntil);
    B_HR19 = B_HR19(:,[n,1:n-1,n+1]);
    Atil = inv(B_HR19);
    A = Atil(1:4,1:4);
    A = A./[A(1,1),A(2,2),A(3,1),A(4,4)]';
    alpha_qp_pr(i) = - A(1,3); 
    probs_HR20(2,i) = alpha_qp_pr(i)<0.04 ;
    probs_BH19(2,i) = student_pos_prior(alpha_qp_pr(i) ,c_alpha_qp,sigma_alpha_qp,nu_alpha_qp);
    %other: 
    alpha_yp(i) = - A(2,3);
    beta_qy(i) = - A(3,2);
    beta_qp(i) = - A(3,3);
    chi(i) = - 1/A(3,4);

end


 

%% Panel A
 
disp('Panel A: Posterior quantiles (0.16,0.5,0.84) of supply elasticity (R1)')
disp(round(quantile(alpha_qp_R1,[0.16,0.5,0.84]),3)) 
disp('Panel A: Posterior quantiles (0.16,0.5,0.84) of supply elasticity (R2)')
disp(round(quantile(alpha_qp_R2,[0.16,0.5,0.84]),3))

%% Panel B
J=10; 
probs_HR20J = squeeze(mean(reshape(probs_HR20,[2,size(probs_HR20,2)/J,J]),2));
probs_BH19J = squeeze(mean(reshape(probs_BH19,[2,size(probs_BH19,2)/J,J]),2));
BFHR20se = std(2*log(probs_HR20J(2,:)./probs_HR20J(1,:)))./sqrt(J);
BFBH19se = std(2*log(probs_BH19J(2,:)./probs_BH19J(1,:)))./sqrt(J); 
ProbsHR20 = mean(probs_HR20,2);
ProbsBH19 = mean(probs_BH19,2); 
BFHR20 = 2*log(ProbsHR20(2)/ProbsHR20(1));
BFBH19 = 2*log(ProbsBH19(2)/ProbsBH19(1));  
disp('Panel B: Bayes factors testing restrictions on alpha(qp): BH19')
disp([round(ProbsBH19,3)',BFBH19,BFBH19se])
disp('Panel B: Bayes factors testing restrictions on alpha(qp): HR20')
disp([round(ProbsHR20,3)',BFHR20,BFHR20se])


%% Panel C: Contribution of supply shocks to the FEVD of the real price of oil
clear; 
load('results/output_BH19_R1.mat') 
load('results/output_BH19_R2_HR20.mat')  
load('results/output_BH19_R2_BH19.mat')    
n = 4; ntil = 5; h = 24; nthin = 10; 
p = 12; variable = 3; shock = 1;
nrep = output_R1.input.nrep; 
FEVDS = zeros(h+1,nrep/nthin,3);
a = 1;
for i = 1:nthin:nrep   
    [~,fevdsR1] = IRF(output_R1.alphas(:,i),output_R1.vecBs(:,i),p,ntil,h,n,0); 
    B_HR19 = reshape(output_BH19_R2_HR20.vecBs(:,i),ntil,ntil); 
    B_BH19 = reshape(output_BH19_R2_BH19.vecBs(:,i),ntil,ntil);   
    [~,fevdsR2_HR] = IRF(output_BH19_R2_HR20.alphas(:,i),vec(B_HR19(:,[n,1:n-1,n+1])),p,ntil,h,n,0);
    [~,fevdsR2_BH] = IRF(output_BH19_R2_BH19.alphas(:,i),vec(B_BH19(:,[n,1:n-1,n+1])),p,ntil,h,n,0);  
    FEVDS(:,a,1) = fevdsR1(:,(shock*n-n)+variable); 
    FEVDS(:,a,2) = fevdsR2_BH(:,(shock*n-n)+variable);
    FEVDS(:,a,3) = fevdsR2_HR(:,(shock*n-n)+variable);
    a = a+1;
end 
hors = [0,24]+1;
disp('Panel C: Contribution of supply shocks to the FEVD of the real price of oil')
disp('Model R1, h=0 and h=24, posterior quantiles (0.16,0.5,0.84)')
disp(round(quantile(FEVDS(hors,:,1),[.16,.5,.84],2),3))
disp('Model R2 + BH19, h=0 and h=24, posterior quantiles (0.16,0.5,0.84)')
disp(round(quantile(FEVDS(hors,:,2),[.16,.5,.84],2),3))
disp('Model R2 + HR20, h=0 and h=24, posterior quantiles (0.16,0.5,0.84)')
disp(round(quantile(FEVDS(hors,:,3),[.16,.5,.84],2),3))
  
