clear; clc;
restoredefaultpath
addpath('data','functions','models','results')
load('dataset_ACR.mat')
%% Model Setup
varnames = {' GDP', ' GDPDEF' , ' CPRINDEX', ' TRARR', ' BOGNONBR', ' FEDFUNDS'};
instrument = {'RR'};
%% Load Dataset
dates = ACR.dates;
Data_Raw  = [100*log(ACR.monthly_GDP), 100*log(ACR.monthly_GDPDEF),...
    100*log(ACR.CPRINDEX), 100*log(ACR.TRARR), 100*log(ACR.BOGNONBR), ACR.FEDFUNDS];
m = ACR.RR;
[~, n] = size(Data_Raw);
idx_nans = sum(isnan([Data_Raw, m]),2)==0;
Data_Raw = Data_Raw(idx_nans,:);
m = m(idx_nans,:); m(isnan(m))=0;
dates = dates(idx_nans);
k = size(m,2);
%% Baseline
input.p = 12; % lags
input.c = 0; % no constant (see ACR code)
input.y = Data_Raw;
input.m = m;
input.exo = 1;
input.lam = 0.5;
input.prior_pers = [ones(1,n),zeros(1,k)];
input.pr_B = @(vB,S_b)priorBadd_ACR( vB,S_b ); % no further prior
input.n_T = 3*12;
input.priorind = 1;
%% Identification restrictions:
ntil = n + k;
for restriction = 1:3
    ZR = NaN(ntil,ntil,2);
    SR.SIGN = NaN(ntil,ntil,2);
    if restriction == 1 % Model R1 (RR as instrument)
        input.nburn = 1000; input.nrep = 5000;
        input.draw_prior = 0;
        shockssel = n;
        ZR(end,1:n-1,2)=1;
        input.ZR = ZR;
    elseif restriction == 2 % Model R2 (Jonas Arias et al. restrictions)
        input.draw_prior = 1; % Here we draw from the prior as we will use it to test IV validity
        shockssel = 1;
        input.nburn = 1000; input.nrep = 50000;
        SR.SIGN(:,shockssel,1) = [ -1, -1, NaN, NaN, NaN,  +1 , NaN(1,k)]'; % Restrictions on A0=inv(B) (contemporanous relations)
        SR.SIGN(:,shockssel,2) = [ NaN, NaN, NaN, NaN, NaN,  +1 , NaN(1,k)]'; % Restrictions on B (impact matrix)
        ZR(4,shockssel,1) = 1; % Not part of systematic component
        ZR(5,shockssel,1) = 1; % Not part of systematic component
        input.SRinstru = @(B)sr_r2(B,n,shockssel); %could bound
    elseif restriction == 3  % Model R3 (combined identification)
        input.nburn = 1000; input.nrep = 5000;
        input.draw_prior = 0;
        shockssel = 1;
        ZR(4,shockssel,1) = 1; % Not part of systematic component
        ZR(5,shockssel,1) = 1; % Not part of systematic component
        SR.SIGN(:,shockssel,1) = [ -1, -1, NaN(1,3),  +1 , NaN(1,k)]'; % Restrictions on A0=inv(B) (contemporanous relations)
        SR.SIGN(:,shockssel,2) = [ NaN, NaN, NaN(1,3),  +1 , NaN(1,k)]'; % Restrictions on B (impact matrix)
        input.SRinstru = @(B)sr_r3(B, k, .5, shockssel); 
    end
    input.ZR = ZR;
    input.targetshock = shockssel;
    [ input.Z1, input.S1, input.Z2, input.S2 ] = get_Z_S( ZR, SR, n, k );
    %% Estimate the Model in the function SVAR_proxyaug
    display(strcat('Estimating Model #',num2str(restriction)'))
    output  = SVAR_BB2021( input );
    save(strcat('results/ACRr',num2str(restriction)),'output')
end


