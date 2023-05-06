clear; clc; 
addpath('functions','results');
load('SRIV_Plain_main')
elas_post = zeros(3,size(output.alphas,2));
elas_pr = zeros(3,size(output.alphas,2));
ntil = 5; n = 4; 
probs_HR20 = zeros(2,size(output.alphas,2));
probs_BH19 = zeros(2,size(output.alphas,2));
probs_HR20b = zeros(2,size(output.alphas,2));
probs_BH19b = zeros(2,size(output.alphas,2));
c_alpha_qp = .1; sigma_alpha_qp = 0.1; nu_alpha_qp = 3;
for irep = 1:size(output.alphas,2)
   % Evaluate at posterior
   B = reshape(output.vecBs(:,irep),ntil,ntil); % posterior draws
   elas_post(1,irep) = B(1,n)/B(3,n);
   elas_post(2,irep) = B(1,1)/B(3,1);
   elas_post(3,irep) = B(1,2)/B(3,2);
   probs_HR20(1,irep) = and(elas_post(2,irep)<0.04,elas_post(3,irep)<0.04); 
   probs_BH19(1,irep) =  exp(prior_add_BH19( output.S_b*vec(reshape(output.vecBs(:,irep),ntil,ntil)), output.S_b  )); 
   % Evaluate at prior
   B = reshape(output.vecBs_pr(:,irep),ntil,ntil); % prior draws
   elas_pr(1,irep) = B(1,n)/B(3,n);
   elas_pr(2,irep) = B(1,1)/B(3,1);
   elas_pr(3,irep) = B(1,2)/B(3,2);
   probs_HR20(2,irep) = and(elas_pr(2,irep)<0.04,elas_pr(3,irep)<0.04);
   probs_BH19(2,irep) =  exp(prior_add_BH19( output.S_b*vec(reshape(output.vecBs_pr(:,irep),ntil,ntil)), output.S_b  )); 
end
   


% Panel A:
inputTable.data = quantile(elas_post(2:3,:),[.16,.5,.84],2);
inputTable.tableRowLabels = {'$\eta_1$','$\eta_2$'};
inputTable.tableColLabels = {'$16\%$',  '$50\%$', '$84\%$'};
inputTable.dataFormat = {'%.2f',3}; % three digits precision for first two columns, one digit for the last
inputTable.tableLabel = 'Panel A: Posterior under R1 and R2';
latexA = latexTable(inputTable); 


% Panel B:
J = 10; 
probs_HR20J = squeeze(mean(reshape(probs_HR20,[2,size(probs_HR20,2)/J,J]),2));
probs_BH19J = squeeze(mean(reshape(probs_BH19,[2,size(probs_BH19,2)/J,J]),2));
BFHR20se = std(2*log(probs_HR20J(2,:)./probs_HR20J(1,:)))./sqrt(J);
BFBH19se = std(2*log(probs_BH19J(2,:)./probs_BH19J(1,:)))./sqrt(J); 
ProbsHR20 = mean(probs_HR20,2);
ProbsBH19 = mean(probs_BH19,2); 
BFHR20 = 2*log(ProbsHR20(2)/ProbsHR20(1));
BFBH19 = 2*log(ProbsBH19(2)/ProbsBH19(1));  
data_B = [[ProbsBH19';ProbsHR20'],[BFBH19;BFHR20],[BFBH19se;BFHR20se]]; 
inputTable.data = data_B;
inputTable.tableRowLabels = {'BH19','HR20'};
inputTable.tableColLabels = {'$\E_{\theta|\tilde{Y}}[p_{2}(\theta)]$', '$ \E_{\theta}[p_{2}(\theta)] $' , '2\ln\widehat{\mbox{BF}}_{10}', 's.e.' };
inputTable.dataFormat = {'%.3f',4}; % three digits precision for first two columns, one digit for the last
inputTable.tableLabel = 'Panel B: Bayes factors testing restrictions on $\eta_1$ and $\eta_2$';
latexB = latexTable(inputTable); 




