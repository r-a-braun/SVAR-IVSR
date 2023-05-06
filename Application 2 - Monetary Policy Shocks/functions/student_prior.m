function p=student_prior(x, c, sigma, nu)

p = (1/sigma)*tpdf((x-c)/sigma,nu);