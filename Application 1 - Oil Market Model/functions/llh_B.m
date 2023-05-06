function llh = llh_B(vB,S,S_b,T,ntil)
% Likelihood for B
B = reshape(S_b'*vB,ntil,ntil);
llh = - T*LogAbsDet(B) -0.5*trace((B\S)/B');
end
