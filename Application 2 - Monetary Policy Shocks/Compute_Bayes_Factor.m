% clear;clc;
load('ACRr2.mat')
% load('output_ACR_rep.mat')
p = output.input.p;
%% Test the IV exclusion restrictions via Bayes Factor
y = output.input.y; % VAR Data of size [T,n]
[T, n]=size(y);
m = output.input.m; % Instruments of size [T,k]
k = size(m,2);
idx_nan = isnan(m); % find NaNs (missings)
m(idx_nan ) = 0; % set NaNs with 0's
Ytilraw = [y, m];
[Tp, ntil] = size(Ytilraw);
T = Tp - p; 
lmX = lagmatrix(Ytilraw,1:p);
Xtil = [ones(T,output.input.c), lmX(p+1:end,:)];
Ytil = Ytilraw(p+1:end,:); 
Xtil = Xtil(output.input.n_T+1:end,:);
Ytil = Ytil(output.input.n_T+1:end,:);  
S0 = output.S0  ;
% S0 =blkdiag(output.S0(1:n,1:n),output.S0(n+1:end,n+1:end));
v0 = output.input.n_T;
Sr = eye(n*k); Sr(1,:)=[];
Sf = eye(n*k); Sf(2:end,:)=[];
M = Ytil(:,n+1)';    
S11 = S0(1:n,1:n); 
S21 = S0(n+1:end,1:n);
S22 = S0(n+1:end,n+1:end); 
logprior_den = zeros(output.input.nrep,1);
logpost_den = zeros(output.input.nrep,1);

for i=1:output.input.nrep 
   Btil = reshape(output.vecBs(:,i),n+1,n+1);
   Sigeta = Btil(n+1:end,n+1:end)*Btil(n+1:end,n+1:end)';
   B = Btil(1:n,1:n);
   Phi = Btil(n+1:end,1:n); 
   Alpha = reshape(output.alphas(:,i),p*ntil+output.input.c,ntil); 
   E = (Btil\(Ytil - Xtil*Alpha)');
   E = E(1:n,:);
   %% Evaluate posterior draws: 
   Spost = S0 + (Ytil - Xtil*Alpha)'*(Ytil - Xtil*Alpha);
   S11p = Spost(1:n,1:n);
   S21p = Spost(n+1:end,1:n);
   S22p = Spost(n+1:end,n+1:end);
   postm =   vec(S21p*inv(S11p)*B) ;  
   postV = kron(B'*inv(S11p)*B,Sigeta); 
   mu_r = Sr*postm + ...
       (Sr*postV*Sf')*inv(Sf*postV*Sf')*(Sf*vec(Phi) - Sf*postm);
   V_r = Sr*postV*Sr' - (Sr*postV*Sf')*inv(Sf*postV*Sf')*(Sf*postV*Sr');   
   logpost_den(i) = logmvnpdf(zeros((n-k)*k,1)',mu_r',V_r);
   
   
   %% Evaluate prior draws: 
   Btil  = reshape(output.vecBs_pr(:,i),n+1,n+1);  
   Sigeta = Btil(n+1:end,n+1:end)*Btil(n+1:end,n+1:end)';
   B = Btil(1:n,1:n);
   Phi = Btil(n+1:end,1:n); 
   postm =  vec( S21*inv(S11)*B ) ;
   postV = kron(B'*inv(S11)*B,Sigeta); 
   mu_r = Sr*postm + ...
       (Sr*postV*Sf')*inv(Sf*postV*Sf')*(Sf*vec(Phi) - Sf*postm);
   V_r = Sr*postV*Sr' - (Sr*postV*Sf')*inv(Sf*postV*Sf')*(Sf*postV*Sr');  
   logprior_den(i) = logmvnpdf(zeros((n-k)*k,1)',mu_r',V_r);
   clc; disp(i/output.input.nrep )
end
% save('results/ACRr2bf','logprior_den','logpost_den')


J = 10; 
maxpost = max(logpost_den);
logpost = log(mean(exp(logpost_den-maxpost))) + maxpost;
logposts = log(mean(reshape(exp(logpost_den-maxpost),J,output.input.nrep/J),2)) + maxpost;
maxpr = max(logprior_den);
logprior = log(mean(exp(logprior_den-maxpr))) + maxpr;
logpriors = log(mean(reshape(exp(logprior_den-maxpr),J,output.input.nrep/J),2)) + maxpr;
disp( 2*(logprior-logpost) )
disp( std(2*(logpriors-logposts))./sqrt(J) )

 
M = output.input.nrep ;
plot(1:M,logprior_den,'b',1:M,logpost_den,'r')
 
 