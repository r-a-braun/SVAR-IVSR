function [ p ] = student_neg_prior( x, c, sigma, nu )
%student_pos_prior calculates height of prior truncated to be negative
% at the point x 
p = (1/sigma)*tpdf((x-c)/sigma,nu)./tcdf(-c/sigma,nu);
end

