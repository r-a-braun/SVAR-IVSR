

function zr = zero_restr(vB,ZA,S_b)
ntil = size(ZA,2);
B = reshape(S_b'*vB,ntil,ntil);
A = inv(B)';
zr = A(ZA==1);
end

