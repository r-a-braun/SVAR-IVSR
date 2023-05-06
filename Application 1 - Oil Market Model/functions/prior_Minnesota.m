% This script constructs the asymmetric conjugate prior given the
% hyperparameter values
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

function [mu_p,Var_p] = prior_Minnesota(n,p,prior_pers,ic,kappa,sig2)
ki = ic + n*p; 
mu_p = zeros(ki,n);
Var_p = zeros(ki,n); 
for var_i = 1:n 
    if p>0
    mu_p(ic+var_i,var_i) = prior_pers(var_i);
    else
    end
Vi = zeros(ki,1); 
for j=1:ki 
    l = ceil((j-ic)/n); % lag length
    idx = mod(j-ic,n);  % variable index
    if idx==0
        idx = n;
    end
    if and(ic==1,j==1) % intercept
        Vi(j) = kappa(3); 
    elseif idx == var_i % own lag
        Vi(j) = (kappa(1)*sig2(var_i))/(l^2*sig2(idx));
    else % lag of other variables)
        Vi(j) = (kappa(2)*sig2(var_i))/(l^2*sig2(idx));
    end     
end
 Var_p(:,var_i) = Vi;  
end

%disp('stop');


end