clear; clc; rng(20)
restoredefaultpath
addpath('data','functions','results') 
load('data_BH19.mat') 
%% Specify model 
start = datetime(1975,2,1);% - calmonths(p); % First datapoint including presample value
yraw = [100*diff(log(BH19.Production)),  100*diff(log(BH19.WIP)),...
    100*diff(log(BH19.RAC)),  BH19.Inventories(2:end)];
date_yraw = BH19.date(2:end);
mraw = BH19.kilian08extended(2:end);
y = yraw(date_yraw>=start,:); m = mraw(date_yraw>=start,:);
date = date_yraw(date_yraw>=start,:); 
p = 12;
k = size(m,2);  
[TpP, n] = size(y);
varnames = {' Prod', ' WIP' , ' RPO', ' INV'}; 
instrument = {'kilian08'};
input.exo = 1;
input.lam = 1;  
input.y = y;  
input.m = m;  
%% Baseline   
input.nburn = 1000;
input.nrep = 10000;   
input.p = p; % lags
input.c = 1; % no constant (see ACR code) 
%% Identification restrictions:
% Ordering: 1) AD, 2) CD, 3) ID 4) AS 
ntil = n + k;
% First k shocks are from the SVAR
ZR = NaN(ntil,ntil,2);    % Place ZERO restrictions on   F(THETA)=[A_{0};L_{0} ]
SR.SIGN = NaN(ntil,ntil,2); % Place SIGN restrictions on F(THETA)=[A_{0};L_{0};\ldots] 
input.prior_pers = [zeros(1,n),zeros(1,k)];
supply_shock = 4;
SR.SIGN(:,supply_shock,1) = [+1, NaN, -1, NaN, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
AD_shock = 1;
SR.SIGN(:,AD_shock,1) = [NaN, +1, +1,NaN, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
ZR([1, 4],AD_shock,1) = 1; % Not part of global production equation
CD_shock = 2;
SR.SIGN(:,CD_shock,1) = [+1, -1, +1, -1, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
ID_shock = 3;
SR.SIGN(:,ID_shock,1) = [+1, NaN, NaN, NaN, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
ZR(end,1:n-1,2)=1;
input.ZR = ZR; 
[ input.Z1, input.S1, input.Z2, input.S2 ] = get_Z_S( ZR, SR, n, k );
%% Estimate the Models 
order = 2; % For all models AS is ordered last now (due to the IV restrictions)
input.n_T = 5*12;
input.draw_prior = 1;   
BH19alpha_qp = 0; % restriction at alpha_pq
input.pr_B = @(vB,S_b)prior_BH19( vB, S_b, BH19alpha_qp, order );
input.SRinstru = @(B)chi_BH19(B, order); % non-linear SR on chi (AS first) 
output_BH19_R2  = SVAR_BB2021( input );
save('results\output_BH19_R2','output_BH19_R2') 

%% Model R2 + BH19
input.draw_prior = 0; BH19alpha_qp = 1; % restriction at alpha_pq
input.pr_B = @(vB,S_b)prior_BH19( vB,S_b, BH19alpha_qp, order );
input.SRinstru = @(B)chi_BH19(B, order); % non-linear SR on chi (AS first) 
output_BH19_R2_BH19  = SVAR_BB2021( input );  
save('results\output_BH19_R2_BH19','output_BH19_R2_BH19')  
 
%% Model R2 + HR20
input.draw_prior = 0; BH19alpha_qp = 0; % BH restriction at alpha_pq 
input.pr_B = @(vB,S_b)priorBadd( vB,S_b,BH19alpha_qp,order );   
input.SRinstru = @(B)HR20(B); % HR restriction at alpha_pq    
output_BH19_R2_HR20 = SVAR_BB2021( input ); 
save('results\output_BH19_R2_HR20','output_BH19_R2_HR20')  

 

