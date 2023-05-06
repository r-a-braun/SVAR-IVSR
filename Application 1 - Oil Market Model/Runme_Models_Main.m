clear; clc;
restoredefaultpath
addpath('data','functions','results') 
load('oilmarket_data.mat')

%% Specify model  
p = 13; 
date = data.date;
y = [ 100*log(data.Prod) , 100*log(data.WIP) ,...
    100*log(data.RAC) , 100*log(data.inventsa) ];
m = data.kilian08extended;
m(isnan(m))=0;   
startdate = datetime(1978,10-2,1); enddate = datetime(2018,11,1);
idx_Y = and(date >=  startdate ,date <= enddate ); 
date = date(idx_Y); m = m(idx_Y); y = y(idx_Y,:);
input.n_T = 5*12; 
k = size(m,2);  
[TpP, n] = size(y); 
varnames = {' dProd', ' WIP' , ' RPO', ' dI'}; 
instrument = {'k08'}; 
input.exo = 1;
input.lam = 1;  
input.y = y;  
input.m = m;  
%% Baseline   
input.nburn = 500;
input.nrep = 5000; % DEFAULT VALUE  
input.p = p; % lags
input.c = 1; % no constant (see ACR code) 
%% Identification restrictions:
ntil = n + k;
% First k shocks are from the SVAR
SR.SIGN = NaN(ntil,ntil,2); % Place SIGN restrictions on F(THETA)=[A_{0};L_{0};\ldots] 
input.prior_pers = [zeros(1,n),zeros(1,k)];
flow_supply_shock = n;
SR.SIGN(:,flow_supply_shock,2) = [-1, -1, +1, NaN, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
flow_demand_shock = 1;
SR.SIGN(:,flow_demand_shock,2) = [+1, +1, +1, NaN, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
storage_demand_shock = 2;
SR.SIGN(:,storage_demand_shock,2) = [+1, -1, +1, +1, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
% Ordering: 1) AD, 2) CD, 3) ID 4) AS
ZR = NaN(ntil,ntil,2);    % Place ZERO restrictions on   F(THETA)=[A_{0};L_{0} ]
input.ZR = ZR;
[ input.Z1, input.S1, input.Z2, input.S2 ] = get_Z_S( ZR, SR, n, k );
%% Estimate the Model in the function SVAR_proxyaug 
Restrictions = {'SR','SRIV'};
Elasticity = {'Plain','HR20','BH19'};
for i = 1:size(Restrictions,2)
    for j = 1:size(Elasticity,2)
        model = strcat(Restrictions(i),'_',Elasticity(j));
        disp(model)
        tic;
        input_ij = input;
        if strcmp(Restrictions{i},'SRIV')
            ZR(end,1:n-1,2)=1;
            input_ij.ZR = ZR;
            lbgamma = 0.025;
            [ input_ij.Z1, input_ij.S1, input_ij.Z2, input_ij.S2 ] = get_Z_S( ZR, SR, n, k );
        elseif strcmp(Restrictions{i},'SR')
            input_ij.acceptall = 1; % Under SR only (and Plain/HR, no ARMH necessary) 
            lbgamma = 0;
        end
        if strcmp(Elasticity{j},'Plain')
            input_ij.SRinstru = @(B)R3_plain(B,lbgamma); % could bound
        elseif strcmp(Elasticity{j},'HR20') 
            input_ij.nburn = 1;
            input_ij.SRinstru = @(B)R3_HR20(B,lbgamma); % could bound 
        elseif strcmp(Elasticity{j},'BH19')
            input_ij.acceptall = 0;
            input_ij.SRinstru = @(B)R3_plain(B,lbgamma); % could bound 
            input_ij.pr_B = @(vB,S_b)prior_add_BH19( vB,S_b ); 
        end  
        
        if and(strcmp(Restrictions{i},'SRIV'),strcmp(Elasticity{j},'Plain'))
            input_ij.nrep = 20000; % Lots of draws for Bayes Factors
            input_ij.draw_prior = 1; 
        else
            input_ij.draw_prior = 0; 
        end 
        output  = SVAR_oilmarket( input_ij );  
        toc;
        save(strcat('results\', model{1}, '_main'),'output')
    end
end

 

