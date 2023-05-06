
function lprior = lnpriorBnew2( vB, S_b, S0, v0 , S_rphi)
% Prior for B
ntil = size(S0,2);
n = size(S_rphi,2);
k = ntil - n;
S_11 = S0(1:n,1:n); S_12 = S0(1:n,n+1:n+k);
S_22 = S0(n+1:n+k,n+1:n+k);
S_22dot1 = S_22 - S_12'*(S_11\S_12);
Btil = reshape(S_b'*vB,ntil,ntil);
B = Btil(1:n,1:n);
Sig_eta = Btil(n+1:end,n+1:end)*Btil(n+1:end,n+1:end)';
vPhi = S_rphi*vec(Btil(n+1:end,1:n));
pm = S_rphi*vec(S_12'*(S_11\B));
pV = S_rphi*kron(B'*(S_11\B),  Sig_eta )*S_rphi';
%pV = S_rphi*kron( inv(S_11) , B'*Sig_eta*B)*S_rphi';
dPhim = vPhi-pm;
lprior1 = - (v0+n)*LogAbsDet(B) - .5*trace( S_11*inv(B*B') );
lprior2 = - (v0+k)/2*LogAbsDet(Sig_eta) - .5*trace( S_22dot1*inv(Sig_eta) );
lprior3 = -0.5*LogAbsDet(pV) - 0.5*(dPhim'*(pV\dPhim));
lprior = lprior1 + lprior2 + lprior3;




end