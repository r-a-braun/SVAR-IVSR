function p=beta_prior(x,alpha_k,beta_k)

%p=(x.^(alpha_k-1).*(1-x).^(beta_k-1))./beta(alpha_k,beta_k);

p=betapdf(x,alpha_k,beta_k);