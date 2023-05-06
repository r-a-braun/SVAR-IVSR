
function [vSw,vBp, notfound,count] = pstar_draw(S, T, ntil, k, ztil, stil, Z1, Z2, S1, S2, Dnn, S_b, W, RInstr)
% Generates a draw from the AR-MH proposal s.t. sign and zero restrictions
n = ntil-k;
flag = 0;
count = 0;
notfound = 0;
while flag==0
    Sig = iwishrnd(S,T); %inv(wish(inv(S),T));
    w = [Draw_w(n,ztil(1:n)); Draw_w(k,ztil(n+1:end))];
    vSw = [ vech(Sig); w];
    vBp = f_map(vSw, ntil, k, ztil,Z1, Z2, Dnn, S_b, W);
    % vSw2 = f_map_inv(vBp, ntil, k, ztil, Z1, Z2, S_b, W);
    % vBp2 = f_map(vSw2, ntil, k, ztil,Z1, Z2, Dnn, S_b, W);
    B = reshape(S_b'*vBp,ntil,ntil);
    flag = 1;
    if nargin>13
        if RInstr(B)==1
        else
            flag = 0;
            continue
        end
    end
    G = [inv(B),B']';
    for j = find(stil~=0)'
        if j<=n
            if sum(S1{j}*G(:,j)>0)==length(S1{j}*G(:,j))
            elseif sum(S1{j}*G(:,j)>0)==0
            else
                flag = 0;
                continue
            end
        else
            if sum(S2{j-k}*G(:,j)>0)==length(S2{j-k}*G(:,j))
            elseif sum(S2{j-k}*G(:,j)>0)==0
            else
                flag = 0;
                continue
            end
        end
    end 
    
    count = count + 1;
    if count>4*1e5
        notfound = 1;
        flag = 1;
    end
end


end