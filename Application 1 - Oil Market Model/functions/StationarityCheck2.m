%function to check the stationarity of draws, no constant
function stat = StationarityCheck2(Beta,n,p)
    stat = 1; 
    if size(Beta,2) == n
        b = Beta';
    else
        b = Beta;
    end
    if size(b,2)==n*p
    else
        b = b(:,2:end);
    end
    LP = [eye(n*(p-1)),zeros(n*(p-1),n)];
    A = [b;LP]; 
    if sum(abs(eig(A))<=1.001) == n*p
        stat = 0;
    end

end