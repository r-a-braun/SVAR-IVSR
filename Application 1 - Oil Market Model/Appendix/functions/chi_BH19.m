function [ flag ] = chi_BH19( B, order )
%CHI_BH19 Summary of this function goes here
%   Detailed explanation goes here
Atil = inv(B);
A = Atil(1:4,1:4);
if order == 1
    % % Ordering: 1)AS, 2) AD, 3) CD, 4) ID 4)
    A = A./[A(1,1),A(2,2),A(3,1),A(4,4)]';
    chi = - 1/A(3,4); 
elseif order == 2
    % % Ordering: 1) AD, 2) CD, 3) ID 4) AS
    A = A./[A(1,2),A(2,1),A(3,4),A(4,1)]'; 
    chi = - 1/A(2,4);
end   
flag = 0;
if chi>1 || chi <0
    return
end
flag = 1;
end

