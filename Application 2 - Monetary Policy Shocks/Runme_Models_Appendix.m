clear; clc;
restoredefaultpath
addpath('data','functions','results')
load('dataset_ACR')
%% Model Setup (CORE)
varnames = {' GDP', ' GDPDEF' , 'CPRINDEX', ' TRARR', ' BOGNONBR', ' FEDFUNDS'}; 
instrument = {'RR'};
%% Load Dataset
dates = ACR.dates; m = ACR.RR;
Data  = [100*log(ACR.monthly_GDP), 100*log(ACR.monthly_GDPDEF),...
    100*log(ACR.CPRINDEX), 100*log(ACR.TRARR), 100*log(ACR.BOGNONBR), ACR.FEDFUNDS];
FinVars = [100*log(ACR.realSP500), ACR.mspread, ACR.BAAspread, ACR.CPspread,ACR.EBP];
finvars_name = {' SP500',' mortgage spreads', ' commercial paper spread',' excess bond premium'};
finvars_name_compact = {'SP','MS', 'CPS','EBP'};
p = 12;

%% Baseline  
input.nburn = 500;
input.nrep = 5000;   
input.p = 12; % lags
input.c = 0; % no constant (see ACR code)   
input.exo = 1;
input.lam = 0.5;
input.pr_B = @(vB,S_b)priorBadd_ACR( vB,S_b ); % no further prior
input.n_T = 3*12;
input.priorind = 1;
input.draw_prior = 0;
shockssel = 1; % index of MP shock 
input.prior_pers = [ones(1,size(Data,2) + 1),zeros(1,1)]; 
k = size(m,2);
models = [1,2,4,5]; %Used in the pictures (BAA not necessary!)
for i = 1:size(FinVars,2)
    disp(strcat('Model: ',num2str(i),'out of', num2str(size(FinVars,2))))
    % Set up raw data
    Data_Raw = [Data(:,1:3),FinVars(:,i),Data(:,4:6)];
    m = ACR.RR; dates = ACR.dates;
    varnames = {' GDP', ' GDPDEF' , ' CPRINDEX', finvars_name{i}, ' TRARR', ' BOGNONBR', ' FEDFUNDS'};  
    [~, n]=size(Data_Raw); 
    % Truncate
    idx_nans = sum(isnan([Data_Raw,m]),2)==0;
    Data_Raw = Data_Raw(idx_nans,:);
    m = m(idx_nans,:);
    m(isnan(m))=0;
    dates = dates(idx_nans);
    disp([dates(1),dates(end)])
    input.dates = dates;
    input.varnames = varnames;
    input.y = Data_Raw; 
    input.m = m; 
    %% Identification restrictions:
    ntil = n + k;
    for restriction = 2:3
        ZR = NaN(ntil,ntil,2);    % Place ZERO restrictions on   F(THETA)=[A_{0};L_{0} ]
        SR.SIGN = NaN(ntil,ntil,2); % Place SIGN restrictions on F(THETA)=[A_{0};L_{0};\ldots]
        if restriction ==2  
            SR.SIGN(:,shockssel,1) = [ NaN(1, n-1),  +1 , NaN(1,k)]'; % Restrictions on B (impact matrix)
            SR.SIGN(:,shockssel,2) = [ -1, -1, NaN(1,4),  +1 , NaN(1,k)]'; % Restrictions on A0=inv(B) (contemporanous relations)
            ZR(5,shockssel,1) = 1; % Not part of systematic component
            ZR(6,shockssel,1) = 1; % Not part of systematic component
            input.SRinstru = @(B)sr_r2(B,n,shockssel); %could bound 
        elseif restriction == 3  
            SR.SIGN(:,shockssel,1) = [ NaN(1, n-1),  +1 , NaN(1,k)]'; % Restrictions on B (impact matrix)
            SR.SIGN(:,shockssel,2) = [ -1, -1, NaN(1,4),  +1 , NaN(1,k)]'; % Restrictions on A0=inv(B) (contemporanous relations)
            ZR(5,shockssel,1) = 1; % Not part of systematic component
            ZR(6,shockssel,1) = 1; % Not part of systematic component
            input.SRinstru = @(B)sr_r3(B, k, .5, shockssel);
        end
        [input.Z1, input.S1, input.Z2,input.S2 ] = get_Z_S( ZR,SR,n,k );
        input.targetshock = shockssel;
        input.ZR = ZR;
        %% Estimate the Model in the function SVAR_proxyaug
        output  = SVAR_BB2021( input );
        save(strcat('results/ACRr',num2str(restriction),finvars_name_compact{i}),'output')
    end

end


