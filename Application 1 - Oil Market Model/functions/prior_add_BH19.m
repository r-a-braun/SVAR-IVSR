function logprior = prior_add_BH19( vB, S_b  )
%PRIOR_ADD_BH19 Summary of this function goes here
%   Detailed explanation goes here
B = reshape(S_b'*vB,5,5);

c_alpha_qp = .1; sigma_alpha_qp = 0.1; nu_alpha_qp = 3;

ela_supply = B(1,1)/B(3,1); 
ela_supply2 = B(1,2)/B(3,2); 

d1=student_pos_prior(ela_supply ,c_alpha_qp,sigma_alpha_qp,nu_alpha_qp);
d2=student_pos_prior(ela_supply2,c_alpha_qp,sigma_alpha_qp,nu_alpha_qp);
logprior = log(d1)+log(d2);
end

