function flag = psi_restrictions2(B,n,shock)
% Checks Sign Restrictions and constraints reliability of instruments
A0 = inv(B)'; 
MPeq = - A0(1:n-1,shock)./(A0(n,shock));
flag = 0; 
for i = 1:2 
    if abs(MPeq(i))>4
        return
    end 
end
flag = 1;

end

