function g=Gradp(f,x0,varargin)
% computes the gradient of f evaluated at x
% uses forward gradients. Adjusts for possible differently scaled x by taking percentage increments
% this function is the equivalent to the gradp function of Gauss
% f should return either a scalar or a column vector
% x0 should be a column vector of parameters
f0=feval(f,x0,varargin{:}); 
[T,c]=size(f0);

if size(x0,2)>size(x0,1)
    x0=x0';
end
k=size(x0,1); % number of parameters wrt which one should compute gradient

h=0.0000001; %some small number

g=zeros(T,k); %will contain the gradient
e=eye(k); 
for j=1:k;
    if x0(j)>1; % if argument is big enough, compute relative number   
        f1=feval(f,(x0.*( ones(k,1) +  e(:,j) *h )),varargin{:});    
        g(:,j)=(f1-f0)/(x0(j)*h);    
    else
        f1=feval(f, x0 +  e(:,j) *h ,varargin{:});    
        g(:,j)=(f1-f0)/h;    
    
    end
    
end
