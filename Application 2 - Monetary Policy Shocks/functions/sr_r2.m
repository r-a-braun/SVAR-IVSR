function flag = sr_r2(B,n,shock)
% Checks Sign Restrictions and constraints reliability of instruments
A0 = inv(B)'; 
MPeq = - A0(1:n,shock)./(A0(n,shock));
flag = 0;  
for i = 1:2 
    if abs(MPeq(i))>4%,abs(MPeq(i))<0.25)
        return
    end  
end
flag = 1; 
end

