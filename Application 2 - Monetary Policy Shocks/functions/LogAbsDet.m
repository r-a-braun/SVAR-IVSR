function x = LogAbsDet(X)
% 
% Computes the log of the absolute value of the determinant of the square matrix
% X.
%
% Uses the LU decomposition
%

[~,U,~]=lu(X);
n=size(U,1);
x=0.0;
for i=1:n
    if U(i,i) == 0.0
        x=-inf;
        return;
    end
    x=x+log(abs(U(i,i)));
end
    