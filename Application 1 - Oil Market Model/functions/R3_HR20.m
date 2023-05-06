function [ flag ] = R3_HR20( B , gamma)
%CHI_BH19 Summary of this function goes here
%   Detailed explanation goes here  
flag = 0; 
ntil = size(B,2); k = 1; n = 4;  
Gamma = B(end-k+1:end,ntil-2*k+1:ntil-k);
CholSigm = B(end-k+1:end,end-k+1:end);
RM = inv(CholSigm*CholSigm')*Gamma*Gamma';   
if min(eig(RM))<gamma
    return
end  
% ela_demand_prod = B(1,n)/B(3,n);
% if ela_demand_prod<-0.8
%     return
% end  

ela_supply = B(1,1)/B(3,1);
if ela_supply > 0.04
    return
end  
ela_supply2 = B(1,2)/B(3,2); 
if ela_supply2 > 0.04
    return
end  

% FEVD0rpo = B(3,:).^2./sum(B(3,:).^2);
% if FEVD0rpo(1) < 0.01 % elasticity is not identified
%     return
% end 
% if FEVD0rpo(2) < 0.01 % elasticity is not identified
%     return
% end 
flag = 1;
end

