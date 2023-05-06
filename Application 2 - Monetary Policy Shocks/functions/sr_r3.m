function flag = sr_r3(B,k,c,shock)
% Checks Sign Restrictions and constraints reliability of instruments
flag = 0;
ntil = size(B,2);
n = ntil - k; 
FEVDh0 =  B.^2 ./sum( B.^2 ,2);
for i = 1:k
    FEVDh0shocks = FEVDh0(n+i,1:ntil-k)./sum(FEVDh0(n+i,1:ntil-k));
    [~,idx]= max(FEVDh0shocks);
    if idx~=1
        return
    end
    if FEVDh0shocks(shock) < c 
        return
    end
end 
A0 = inv(B)'; 
MPeq = - A0(1:n,shock)./(A0(n,shock));
for i = 1:2
    if abs(MPeq(i))>4
        return
    end  
end
flag = 1; 
 
end
