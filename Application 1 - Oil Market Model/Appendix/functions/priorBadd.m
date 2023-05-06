function log_priors = priorBadd( vB, S_b, elasupply, order )
%PRIORBADD Summary of this function goes here
%   Detailed explanation goes here
B = reshape(S_b'*vB,5,5);
Atil = inv(B);
A = Atil(1:4,1:4);

% % % Ordering: 1) AS 2) AD, 3) CD, 4) ID 
if order == 1 % Supply shock is first
    A = A./[A(1,1),A(2,2),A(3,1),A(4,4)]';
    
    
    alpha_qp = - A(1,3);
    alpha_yp = - A(2,3);
    beta_qy = - A(3,2);
    beta_qp = - A(3,3);
    chi = - 1/A(3,4);
elseif order == 2
    % Supply shock is last
    % % % Ordering: 1) AD, 2) CD, 3) ID 4) AS 
    A = A./[A(1,2),A(2,1),A(3,4),A(4,4)]';
    alpha_yp = - A(1,3);
    beta_qy = - A(2,2);
    beta_qp = - A(2,3);
    chi = - 1/A(2,4);
    alpha_qp = - A(4,3);
end


% SPECIFY PRIOR PARAMETERS
% alpha(qp): short-run price elasticity of oil supply (sign: positive)
c_alpha_qp = .1; sigma_alpha_qp = 0.1; nu_alpha_qp = 3;
% alpha(yp): short-run oil price elasticity of global demand (sign: negative)
c_alpha_yp = -0.05; sigma_alpha_yp = 0.1; nu_alpha_yp = 3;
% beta(qy): income elasticity of oil demand (sign: positive)
c_beta_qy = 0.7; sigma_beta_qy = 0.2; nu_beta_qy = 3;
% beta(qp): short-run price elasticity of oil demand (sign: negative)
c_beta_qp = -0.1; sigma_beta_qp = 0.2; nu_beta_qp = 3; 
% prior chi: OECD fraction of true oil inventories (about 60-65%)
alpha_k = 15; beta_k = 10;  
% EVALUATE DENSITIES
prior2 = student_neg_prior(alpha_yp,c_alpha_yp,sigma_alpha_yp,nu_alpha_yp);
prior3 = student_pos_prior(beta_qy,c_beta_qy,sigma_beta_qy,nu_beta_qy);
prior4 = student_neg_prior(beta_qp,c_beta_qp,sigma_beta_qp,nu_beta_qp); 
prior5 = beta_prior(chi,alpha_k,beta_k);
if elasupply == 1
    prior1 = student_pos_prior(alpha_qp,c_alpha_qp,sigma_alpha_qp,nu_alpha_qp);
    log_priors =  log(prior1) + log(prior2) + log(prior3) + log(prior4) + log(prior5);
else
    log_priors =  log(prior2) + log(prior3) + log(prior4) + log(prior5);
    
end

end

