function [ flag ] = HR20( B )
%CHI_BH19 Summary of this function goes here
%   Detailed explanation goes here 
Atil = inv(B);
A = Atil(1:4,1:4); 
% % Ordering: 1) AD, 2) CD, 3) ID 4) AS
A = A./[A(1,2),A(2,1),A(3,4),A(4,1)]';
chi = - 1/A(2,4); 
flag = 0;
if chi<0 || chi>1  || -A(end,3)>0.04  || -A(end,3)< 0 
    return
end
flag = 1;
end

