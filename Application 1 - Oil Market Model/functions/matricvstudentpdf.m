
% SUBFUNCTIONS
function logpdf = matricvstudentpdf(X,M,P,Q,v)
% res = matricvstudentpdf(X,M,P,Q,v)
% PURPOSE: evaluate LOG Matricvariate Student density at point X
% reference: Bauwens, Lubrano, Richard (1999) A.2.7
% DEPENDS: multgammaln
[p,q] = size(M);
logCMt = 0.5*p*q*log(pi) - 0.5*v*log(det(Q)) - 0.5*q*log(det(P)) + ...
    multgammaln(q, 0.5*v) - multgammaln(q, 0.5*(v+p));
logkernel = -0.5*(v+p)*log(det(Q + (X-M)'*P*(X-M)));
logpdf = -logCMt + logkernel;
end


function res = multgammaln(N,a)
% PURPOSE: evaluate log of the N-dimensional multivariate Gamma function,
% defined: pi^{N(N-1)/4} \prod_{n=1}^N \Gamma ((2a+1-n)/2) 
nlst = 1:N ;
v1n2 = (2*a+1-nlst)/2 ;
gamln = gammaln(v1n2) ;
sumgamln = sum(gamln) ;
res = 0.25*N*(N-1)*log(pi) + sumgamln;
end
