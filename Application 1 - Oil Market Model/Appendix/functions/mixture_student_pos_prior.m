function [ p ] = mixture_student_pos_prior( x )
%student_pos_prior calculates height of prior truncated to be positive
% at the point x 
theta = 0.8;
sigma = 0.2;
c = 0.1; 
nu = 3;
ub = 0.04; 
if and(x>0,x <=ub)
p = theta/ub + (1-theta)*(1/sigma)*tpdf((x-c)/sigma,nu)./(1 - tcdf(-c/sigma,nu));
elseif and(x>0,x>ub)
p = (1-theta)*(1/sigma)*tpdf((x-c)/sigma,nu)./(1 - tcdf(-c/sigma,nu));
else
   p = 0; 
end
end

