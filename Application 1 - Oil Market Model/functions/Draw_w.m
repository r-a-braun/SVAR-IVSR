
function w = Draw_w(n,ztil)
w = zeros(sum(1:n)-sum(ztil),1);
idx = 0 ;
for j = 1:n
    s = n + 1 - j - ztil(j);
    x1j = randn(s,1);
    w(idx+1:idx+s) = (x1j./norm(x1j));
    idx = idx + s;
end
end



