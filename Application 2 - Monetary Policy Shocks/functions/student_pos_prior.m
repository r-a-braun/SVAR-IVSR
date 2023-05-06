function [ p ] = student_pos_prior( x, c, sigma, nu )
%student_pos_prior calculates height of prior truncated to be positive
% at the point x 
p = (1/sigma)*tpdf((x-c)/sigma,nu)./(1 - tcdf(-c/sigma,nu));
end

