
function vecB = f_map(x, ntil, k, ztil, Z1, Z2, Dnn, S_b, W)
% maps [vech(Sigma), w] to S_B*vec(B)]
n = ntil-k;
nS = size(Dnn,2);
Sigma = reshape(Dnn*x(1:nS),ntil,ntil);
P = chol(Sigma)';
G = [inv(P),P']';
w = x(nS+1:end);
Q1 = zeros(n,n);
idx = 0;
for j = 1:n
    s = n + 1 - j - ztil(j) ;
    wj = w(idx+1:idx+s);
    Mj = [Q1(:,1:j-1)'; Z1{j}*G(:,1:n); W{j,1}];
    [K,R]=qr(Mj');
    for i = n-s+1:n
        if R(i,i)<0
            K(:,i)= -K(:,i);
        end
    end
    Kj = K(:,n-s+1:n);
    Q1(:,j) = Kj*wj;
    idx = idx + s;
end
Q2 = zeros(k,k);
for j = 1:k
    s = k + 1 - j - ztil(n+j);
    Mj = [Q2(:,1:j-1)'; Z2{j}*G(:,n+1:end); W{j,2}];
    wj = w(idx+1:idx+s);
    [K,R]=qr(Mj');
    for i = k-s+1:k
        if R(i,i)<0
            K(:,i)=-K(:,i);
        end
    end
    Kj = K(:,k-s+1:k);
    Q2(:,j) = Kj*wj;
    idx = idx + s;
end
Q = eye(ntil);
Q(1:n,1:n) = Q1;
Q(n+1:end,n+1:end) = Q2;
B = P*Q;
vecB = S_b*vec(B);
end