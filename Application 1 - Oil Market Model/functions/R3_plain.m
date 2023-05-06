function [ flag ] = R3_plain( B,gamma )
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
flag = 1;
end

