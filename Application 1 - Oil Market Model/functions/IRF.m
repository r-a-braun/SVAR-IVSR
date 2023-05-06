function [irfs,fevds] = IRF(valpha,vB,p,ntil,h,n,normalization)
%IRF Summary of this function goes here
%   Detailed explanation goes here
Phi = reshape(valpha,p*ntil+1,ntil);
PhiAR = Phi(2:end,:)';
B =reshape( vB,ntil,ntil);
if normalization==0
B = B.*(-sign(B(1,:)));
elseif normalization==1
B = B.*(sign(B(3,:)));
end
% IRFs
[AA, J] = Companion( PhiAR , 0);
Apoweri = AA;
Phis = zeros(ntil,ntil,h+1);
Phis(:,:,1) = eye(ntil);
for ii = 1:h
    Phis(:,:,ii+1) = J*Apoweri*J';
    Apoweri = Apoweri * AA;
end
[AA, J] = Companion( PhiAR , 0);
Apoweri = AA;
Phis = zeros(ntil,ntil,h+1);
Phis(:,:,1) = eye(ntil);
for ii = 1:h
    Phis(:,:,ii+1) = J*Apoweri*J';
    Apoweri = Apoweri * AA;
end
Theta = zeros(n,n,h+1);
FEVD = zeros(n,n,h+1);
for ii = 0:h
    Theta_ii = Phis(:,:,ii+1)*B;
    Theta(:,:,ii+1) = Theta_ii(1:n,1:n);
    FEVD(:,:,ii+1) = sum(Theta(:,:,1:ii+1).^2,3)./sum(sum(Theta(:,:,1:ii+1).^2,3),2);
end
fevds = reshape(FEVD,n^2,h+1)'; 
irfs = reshape(Theta,n^2,h+1)';  

end

