clear; clc;
n = 6; quants = [.05,.5,.95];

addpath('data','functions','results')
load('ACRr1')
outputR1 = output;
for i = 1:size(outputR1.alphas,2) 
    Btil = reshape(outputR1.vecBs(:,i), n+1, n+1);
    Btil=Btil(:,[n,1:n-1,n+1]); 
    A0 = inv(Btil(1:n,1:n))'; 
    MPeqR1(:,i) = - A0(1:n-1,1)./(A0(n,1)); 
end
TaylorRule(:,:,1) =   quantile(MPeqR1,quants,2);
load('ACRr2') 
outputR2 = output;
for i = 1:size(outputR2.alphas,2) 
    Btil = reshape(outputR2.vecBs(:,i), n+1, n+1); 
    A0 = inv(Btil(1:n,1:n))'; 
    MPeqR2(:,i) = - A0(1:n-1,1)./(A0(n,1)); 
end
TaylorRule(:,:,2) =   quantile(MPeqR2,quants,2);

load('ACRr3')
outputR3 = output;
for i = 1:size(outputR3.alphas,2) 
    Btil = reshape(outputR3.vecBs(:,i), n+1, n+1); 
    A0 = inv(Btil(1:n,1:n))'; 
    MPeqR3(:,i) = - A0(1:n-1,1)./(A0(n,1)); 
end
TaylorRule(:,:,3) =   quantile(MPeqR3,quants,2);


%% COMPUTE IRFS AND POSTERIOR OF TAYLOR RULE COEFFICIENTS
Data = [TaylorRule(1:3,:,1);TaylorRule(1:3,:,2);TaylorRule(1:3,:,3)]; 
%Print Results  
inputTable.data = Data';
inputTable.tableRowLabels = {'$5\%$',  '$50\%$', '$95\%$'};
inputTable.tableColLabels = {'$\xi_y$','$\xi_{\pi}$' , '$\xi_{c}$','$\xi_y$','$\xi_{\pi}$' , '$\xi_{c}$','$\xi_y$','$\xi_{\pi}$' , '$\xi_{c}$'};
inputTable.dataFormat = {'%.2f',9}; % three digits precision for first two columns, one digit for the last
inputTable.tableCaption = 'Posterior policy rule parameters';
inputTable.tableLabel = 'MPrule';
latex = latexTable(inputTable); 
