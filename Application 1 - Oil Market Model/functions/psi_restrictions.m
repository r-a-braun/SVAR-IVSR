function flag = psi_restrictions(B,n,shock)
% Checks Sign Restrictions and constraints reliability of instruments
A0 = inv(B)'; 
MPeq = - A0(1:n-1,shock)./(A0(n,shock));
flag = 0; 
% Bimp = B(:,shock)./B(n,shock);
for i = 1:2 
    if MPeq(i)>4
        return
    end 
%     if or(abs(Bimp(1))>5,abs(Bimp(2))>5)
%         return
%     end
end
flag = 1;

end

