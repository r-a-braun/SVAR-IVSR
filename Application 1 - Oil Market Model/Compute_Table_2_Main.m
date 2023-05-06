clear; clc; addpath('functions','results')
Restrictions = {'SR','SRIV'};
Elasticity = {'HR20','BH19'}; 
ntil = 5; n = 4; hor = 24;
quants = [.16,.5,.84];
IRFS = []; FEVDS = []; Elas = [];
for i = 1:size(Restrictions,2)
    for j = 1:size(Elasticity,2)
        model = strcat(Restrictions(i),'_',Elasticity(j));
        disp(model) 
        load(strcat('results\', model{1}, '_main'),'output')
        irf_draws = zeros(hor+1,n^2,size(output.alphas,2));
        fevd_draws = zeros(hor+1,n^2,size(output.alphas,2));
        elasticities = zeros(3,size(output.alphas,2));
        for irep = 1:size(output.alphas,2)
            B = reshape(output.vecBs(:,irep),ntil,ntil);
            elasticities(1,irep) = B(1,n)/B(3,n);
            elasticities(2,irep) = B(1,1)/B(3,1);
            elasticities(3,irep) = B(1,2)/B(3,2);  
            B = B(:,[n,1:n-1,n+1]);
            [irf_draws(:,:,irep),fevd_draws(:,:,irep)] = IRF(output.alphas(:,irep),vec(B),...
                output.input.p,ntil,hor,n,1);
        end 
        eval(strcat('Elas.',model{1},'=quantile(elasticities,quants,2);'))
        eval(strcat('IRFS.',model{1},'=quantile(irf_draws,quants,3);'))
        eval(strcat('FEVDS.',model{1},'=quantile(fevd_draws,quants,3);')) 
    end
    
end

variable = 3;  
horizons = 1+12;
FEVDsSR = [];
for j = 1:size(Elasticity,2)
    fevdj = [];
    for shock = 1:3 
        eval(strcat('fevdj = [fevdj,squeeze(FEVDS.SR_',Elasticity{j},...
            '(horizons,(shock*n-n)+variable,:))];'))
    end
    FEVDsSR = [FEVDsSR;vec(fevdj)'];
end
FEVDsIVSR = [];
for j = 1:size(Elasticity,2)
    fevdj = [];
    for shock = 1:3 
        eval(strcat('fevdj = [fevdj,squeeze(FEVDS.SRIV_',Elasticity{j},...
            '(horizons,(shock*n-n)+variable,:))];'))
    end
    FEVDsIVSR = [FEVDsIVSR;vec(fevdj)'];
end
 
%% Table  


variable = 3;  
horizons = [1, 1+24];
FEVDsSR = [];
for j = 1:size(Elasticity,2)
    fevdj = [];
    for shock = 1:3 
        eval(strcat('fevdj = [fevdj,squeeze(FEVDS.SR_',Elasticity{j},...
            '(horizons,(shock*n-n)+variable,2))];'))
    end
    fevdj16 = [];
    for shock = 1:3 
        eval(strcat('fevdj16 = [fevdj16,squeeze(FEVDS.SR_',Elasticity{j},...
            '(horizons,(shock*n-n)+variable,1))];'))
    end
    fevdj84 = [];
    for shock = 1:3 
        eval(strcat('fevdj84 = [fevdj84,squeeze(FEVDS.SR_',Elasticity{j},...
            '(horizons,(shock*n-n)+variable,3))];'))
    end
    FEVDsSR = [FEVDsSR;vec(fevdj)';vec(fevdj16)';vec(fevdj84)'];
end


FEVDsIVSR = [];
for j = 1:size(Elasticity,2)
    fevdj = [];
    for shock = 1:3 
        eval(strcat('fevdj = [fevdj,squeeze(FEVDS.SRIV_',Elasticity{j},...
            '(horizons,(shock*n-n)+variable,2))];'))
    end
    fevdj16 = [];
    for shock = 1:3 
        eval(strcat('fevdj16 = [fevdj16,squeeze(FEVDS.SRIV_',Elasticity{j},...
            '(horizons,(shock*n-n)+variable,1))];'))
    end
    fevdj84 = [];
    for shock = 1:3 
        eval(strcat('fevdj84 = [fevdj84,squeeze(FEVDS.SRIV_',Elasticity{j},...
            '(horizons,(shock*n-n)+variable,3))];'))
    end
    FEVDsIVSR = [FEVDsIVSR; vec(fevdj)';vec(fevdj16)';vec(fevdj84)'];
end



 
Data = [FEVDsSR; FEVDsIVSR]; 
inputTable=[];
%Print Results  
inputTable.data = Data;
inputTable.tableColLabels= {'$h=0$',  '$h=24$','$h=0$', '$h=24$','$h=0$', '$h=24$'};
inputTable.tableRowLabels  = {'$SR+HR20$' , '$SR+HR20 (16)$' , '$SR+HR20 (84)$' ,...
    '$SR+BH19$' ,'$SR+BH19(16)$','$SR+BH19(84)$',...
    '$SRIV+HR20$','$SRIV+HR20(16)$','$SRIV+HR20(84)$', ...
    '$SRIV+BH19$','$SRIV+BH19(16)$','$SRIV+BH19(84)$'};
inputTable.dataFormat = {'%.2f',6}; % three digits precision for first two columns, one digit for the last
inputTable.tableCaption = 'Forecast Error Variance Decomposition (h=16)';
inputTable.tableLabel = 'FEVD';
latex = latexTable(inputTable); 

 

