function v = vech(x) 
v = x(tril(ones(size(x,2)))==1); 
