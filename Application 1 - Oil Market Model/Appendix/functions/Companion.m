function [ A, J , nu] = Companion( Bhat , inc)
% Input: [c, A_1,...,A_p] of size (K, N*K+1)
% Output: sparse A, VAR(1) form
[K, Kpinc]=size(Bhat);
p = (Kpinc-inc)/K; 
nu = zeros(K*p,1);
if inc == 1
    nu(1:K) = Bhat(:,1);
end  
A = [Bhat(:,1+inc:end);[speye(K*(p-1)),sparse(K*(p-1),K)]]; 
J = [speye(K),sparse(K,K*(p-1))];   
