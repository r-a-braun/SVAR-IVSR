
function vPcholw = f_map_inv(vB, ntil, k, ztil, Z1, Z2, S_b, W)
%vecB = f_g(x, ntil, k, ztil, Z1, Z2, Dnn, S_b, W)
%F_H Summary of this function goes here
%   Detailed explanation goes here
n = ntil-k;
B = reshape(S_b'*vB,ntil,ntil);
[Qp,R]=qr(B');
for i = 1:size(R,2)
    if R(i,i)<0
        Qp(:,i) = - Qp(:,i);
    end
end
Pchol = B*Qp;
Q1 = Qp(1:n,1:n)';
Q2 = Qp(n+1:end,n+1:end)';
G = [inv(Pchol),Pchol']';
w = [];
idx = 0;
for j = 1:n
    s = n + 1 - j - ztil(j) ;
    Mj = [Q1(:,1:j-1)'; Z1{j}*G(:,1:n); W{j,1}];
    [K,R]=qr(Mj');
    for i = n-s+1:n
        if R(i,i)<0
            K(:,i)= -K(:,i);
        end
    end
    Kj = K(:,n-s+1:n);
    wj = Kj\Q1(:,j);
    w = [w;wj];
    idx = idx + s;
end
for j = 1:k
    s = k + 1 - j - ztil(n+j);
    Mj = [Q2(:,1:j-1)'; Z2{j}*G(:,n+1:end); W{j,2}];
    [K,R]=qr(Mj');
    for i = k-s+1:k
        if R(i,i)<0
            K(:,i)=-K(:,i);
        end
    end
    Kj = K(:,k-s+1:k);
    wj = Kj\Q2(:,j);
    w = [w; wj];
    idx = idx + s;
end


vPcholw = [ vech(Pchol*Pchol');  w];
end