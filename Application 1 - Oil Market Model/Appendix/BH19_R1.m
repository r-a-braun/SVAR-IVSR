clear; clc;
restoredefaultpath
addpath('data','functions','results') 
load('data_BH19.mat') 
%% Specify model
p = 12;
start = datetime(1975,2,1);% - calmonths(p); % First datapoint including presample value
yraw = [100*diff(log(BH19.Production)),  100*diff(log(BH19.WIP)),...
    100*diff(log(BH19.RAC)),  BH19.Inventories(2:end)];
date_yraw = BH19.date(2:end);
m = BH19.kilian08extended;
y = yraw(date_yraw>=start,:); m = m(date_yraw>=start,:);
date = date_yraw(date_yraw>=start,:);
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
%% Estiamte the Models  
ntil = n + k;
ZR = NaN(ntil,ntil,2);    % Place ZERO restrictions on   F(THETA)=[A_{0};L_{0} ]
SR.SIGN = NaN(ntil,ntil,2); % Place SIGN restrictions on F(THETA)=[A_{0};L_{0};\ldots] 
supply_shock = 1;
SR.SIGN(:,supply_shock,1) = [+1, NaN, -1, NaN, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
ZR([2, 4],supply_shock,1) = 1; % Not part of supply equation
AD_shock = 2;
SR.SIGN(:,AD_shock,1) = [NaN, +1, +1,NaN, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
ZR([1, 4],AD_shock,1) = 1; % Not part of global production equation
CD_shock = 3;
SR.SIGN(:,CD_shock,1) = [+1, -1, +1, -1, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
ID_shock = 4;
SR.SIGN(:,ID_shock,1) = [+1, NaN, NaN, NaN, NaN]'; % Restrictions on A0=inv(B) (contemporanous relations)
% Ordering: 1) AS, AD, 2) CD, 3) ID  
input.ZR = ZR; 
[ input.Z1, input.S1, input.Z2, input.S2 ] = get_Z_S( ZR, SR, n, k );
input.pr_B = @(vB,S_b)prior_BH19( vB, S_b, 1, 1 ); % elasupply gets a prior, order=1 (AS first)
input.SRinstru = @(B)chi_BH19(B, 1); % non-linear SR on chi (AS first)
input.draw_prior = 0; 
input.n_T = 5*12; 
input.prior_pers = [zeros(1,n),zeros(1,k)];
output_R1  = SVAR_BB2021( input ); 
save('results\output_BH19_R1','output_R1')

% % 
% % 
% for i = 1:output_R1.input.nrep
%     B = reshape(output_R1.vecBs(:,i),5,5);
%     Atil = inv(B);
%     A = Atil(1:4,1:4);
%     A = A./[A(1,1),A(2,2),A(3,1),A(4,4)]';
%     alpha_qp(i) = - A(1,3);
%     alpha_yp(i) = - A(2,3);
%     beta_qy(i) = - A(3,2);
%     beta_qp(i) = - A(3,3);
%     chi(i) = - 1/A(3,4);  
% end
% 

