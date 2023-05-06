function output = proxy_augmented_SVAR( input )
%SVAR_PROXYAUG Summary of this function goes here
%   Detailed explanation goes here
%% Load Input:
y = input.y; % VAR Data of size [T,n]
m = input.m; % Instruments of size [T,k]
idx_nan = isnan(m); % find NaNs (missings)
m(idx_nan ) = 0; % set NaNs with 0's
p = input.p; % VAR lag order
n = size(y,2);
k = size(m,2);
Ytilraw = [y, m];
[Tp, ntil] = size(Ytilraw);
T = Tp - p;
output = [];

%% Set Regressor Matrix
lmX = lagmatrix(Ytilraw,1:p);
if input.c == 0
    Xtil = lmX(p+1:end,:);
else
    Xtil = [ones(T,1),lmX(p+1:end,:)];
end
Ytil = Ytilraw(p+1:end,:);
vYtil = vec(Ytil);

%% Set zero Constraints on vec(A)
Alre = zeros(ntil);
Alre(1:n,1:n) = ones(n);
if input.exo == 1
    if input.c == 0
        RestrLag = repmat(Alre,p,1) ;
    else
        RestrLag = [ones(1,ntil);repmat(Alre,p,1)];
    end
else
    RestrLag = [ones(input.c,ntil);repmat(ones(ntil),p,1)];
end
[na1,na2]=size(RestrLag);
Ika = speye(na1*na2);
Ika(vec(RestrLag)==0,:)=[];

%% Set identifying constraints on vec(B) and other functions
B_r = ones(ntil);
B_r(input.ZR(:,:,2)==1)  = 0;
S_b = speye(ntil^2);
S_rphi = eye(n*k);
S_rphi(vec(B_r(n+1:end,1:n))==0,:)=[];

S_b(vec(B_r==0),:)=[];
ztil = zeros(ntil,1);
W = cell(ntil,2); % gets random W
stil = zeros(ntil,1);
for j=1:n
    ztil(j) = size(input.Z1{j},1);
    stil(j) = size(input.S1{j},1);
    W{j,1} = randn(n-(j-1+ztil(j)),n);
end
for j=1:k
    ztil(j+n) = size(input.Z2{j},1);
    stil(j+n) = size(input.S2{j},1);
    W{j,2} = randn(k-(j-1+ztil(j+n)),k);
end
Dnn = duplication(ntil);  % Duplication   matrix
if isfield(input,'SRinstru')==1 % Training sample exists:
    RInstr = input.SRinstru;
end
%% Set Prior distributions:
nA = size(Ika,1);
if isfield(input,'ytrain')==1 % Training sample exists:
    if input.c == 0
        Xtr = lagmatrix(input.ytrain,1:p);
    else
        Xtr = [ones(length(input.ytrain),1),lagmatrix(input.ytrain,1:p)];
    end
    Xtr = Xtr(p+1:end,:);
    Ytr = input.ytrain(p+1:end,:);
    Ztr = kron(eye(ntil),Xtr)*Ika';
    Atrain = reshape(Ika'*((Ztr'*Ztr)\(Ztr'*vec(Ytr))),ntil*p+input.c,ntil);
    U = Ytr - Xtr*Atrain;
    Str = U'*U;
    vtr = length(input.ytrain) - ntil*p-1 ;
    V_0 = input.inflator*(Ika*inv(kron(inv(Str/vtr),Xtr'*Xtr))*Ika'); % Just Weakly informative
    iV_0 = inv(V_0);
    v0 = vtr/input.inflator;
    a_0 = Ika*vec(Atrain);
    iV0a0 = iV_0*a_0;
    S0 = (Str/vtr)*v0;
else % choose some uninformative values
    %Sig_phi = S_rphi*kron(eye(k),input.Vphi)*S_rphi'; 
    S0 = eye(n); v0 = n+2;
    S_eta = eye(k); 
    v_eta = k+2;
    %S0 = eye(ntil); v0 = ntil+1; 
    a_0 = zeros(nA,1); V_0 = 100*eye(nA); 
    sig2vars = zeros(n+k,1);
    for i = 1:n+k
        Yii=Ytil(2:end,i); Xii = [ones(size(Yii,1),1),Ytil(1:end-1,i)];
        bii = (Xii'*Xii)\(Xii'*Yii);
        sig2vars(i)= mean( (Yii-Xii*bii).^2 );
    end
    kappa = [input.lam,input.lam/3,1e04] ;
    prior_pers = [ones(1,n),zeros(1,k)];
    [mu_p,Var_p] = prior_Minnesota(n+k,p,prior_pers,input.c,kappa,sig2vars);
    
    a_0 = Ika*vec(mu_p);
    V_0 = sparse(Ika*diag(vec(Var_p))*Ika'); 
    iV_0 = inv(V_0); iV0a0 = iV_0*a_0;
     
    S0aug = zeros(ntil); 
    S0aug(1:n,1:n) = S0; 
    S0aug(n+1:end,n+1:end) = S_eta;
    vaug = ntil+2;
    kappa_phi = 1e05;%input.kappa_phi;
end
% lpriorB = @(vB)lnpriorB( vB, S_b, S0aug, vaug ); % + input.logpvB2(vB)
%lpriorB = @(vB)lnpriorBnew( vB, S_b, S0, v0, kappa_phi, S_eta, v_eta,S_rphi); % + input.logpvB2(vB)
lpriorB = @(vB)lnpriorBnew2( vB, S_b, S0aug, vaug , S_rphi); % + input.logpvB2(vB)

%% Define some functions needed during the algorithm
v_prop = T + vaug ;
if isfield(input,'SRinstru')==1 % Training sample exists:
    proposal_draw = @(xS)pstar_draw(xS, v_prop, ntil, k, ztil, stil, input.Z1, input.Z2, input.S1, input.S2, Dnn, S_b, W, RInstr);
else
    proposal_draw = @(xS)pstar_draw(xS, v_prop, ntil, k, ztil, stil, input.Z1, input.Z2, input.S1, input.S2, Dnn, S_b, W);
end
logtarget = @(xvecB, xSSE) llh_B(xvecB,xSSE,S_b,T,ntil)  + lpriorB(xvecB) ;  % not normalized
finv = @(xvB)f_map_inv(xvB, ntil, k, ztil, input.Z1, input.Z2, S_b, W); % Mapping from vec(B) to [vech(Sigma)',w']'
logproposal = @(xvPw ,xvB , xS )pstar_eval(xvPw, xvB, Dnn, xS, v_prop, finv, input.ZR(:,:,1), S_b) ;

%% Initialize the MCMC algorithm
Ztil = kron(speye(ntil),Xtil)*Ika';
alpha_ols  = (Ztil'*Ztil)\(Ztil'*vYtil);
Alpha_ols = reshape(Ika'*alpha_ols,na1,ntil); % OLS
Sig_ols = (Ytil-Xtil*Alpha_ols)'*(Ytil-Xtil*Alpha_ols)/T; % OLS
iV_p = iV_0 + Ika*kron(inv(Sig_ols),Xtil'*Xtil)*Ika';
alpha_hat = iV_p\(iV0a0 + Ika*vec(Xtil'*Ytil*inv(Sig_ols)));  % Posterior type of mean
alpha = alpha_hat;
Alpha = reshape(Ika'*alpha,na1,ntil);
U = Ytil - Xtil*Alpha;
S = U'*U ;
S_prop = S + S0aug;
ninit = 500;
IS_w = zeros(1,ninit);
vBpcands = zeros(size(S_b,1),ninit);
for j = 1:ninit % Choose a high density point of this density
    [vSwp, vBp] = proposal_draw(S_prop);
    IS_w(j) =   logtarget(vBp,S) -  logproposal(vSwp,vBp,S_prop) ;
    vBpcands(:,j) =  vBp; 
end
var(IS_w)
[~,idx_max] =  max(IS_w);
vB = vBpcands(:,idx_max);
vSw = finv(vB);
B = reshape(S_b'*vB,ntil,ntil); 
iSig = inv(B*B');
Sig = B*B';

%% Save the following quantities:
alpha_draws = zeros(na1*na2,input.nrep);
vecB_draws = zeros(ntil^2,input.nrep);
%% MCMC Algorithm:
logc  = 0;
acc_rate = 1;
it_print = 100;
for irep = 1:input.nrep + input.nburn
    if mod(irep,it_print) == 0
        clc;
        disp([irep/(input.nrep + input.nburn),acc_rate/irep]);
    end
    
    %% Draw Alpha | Sigma, Y 
    iV_p = iV_0 + Ika*kron(iSig,Xtil'*Xtil)*Ika';
    p_m = iV_p\(iV0a0 + Ika*vec(Xtil'*Ytil*iSig));
    stat = 1; countstat = 1;
    while and(stat==1,countstat<1000)
        alphac = p_m + chol(iV_p,'lower')'\randn(nA,1);
        Alphac = reshape(Ika'*alphac,na1,ntil);
        if p>0
            stat = StationarityCheck2(Alphac,ntil,p);
        else
            stat = 0;
        end
        if stat == 0
            Alpha = Alphac; 
            alpha = alphac;
        end
        countstat = countstat + 1;
    end   
    S = (Ytil - Xtil*Alpha)'*(Ytil - Xtil*Alpha);
    
    %% Draw B | alpha, Y  by Accept Reject - Metropolis Hastings (AR-MH) step
    S_prop = S + S0aug;  % Update Proposal distribution for B
    if irep<input.nburn % Tune the constant during the burn-in...
        logc = (irep-1)/irep*logc + 1/irep*( logtarget(vB,S) - logproposal(vSw, vB,  S_prop )) ;
    end
    flag = 0;
    count0 = 0;
    while flag == 0
        count0 = count0 + 1;
        [vSwp, vBp, notfound ] = proposal_draw(S_prop);
        alpARc = logtarget(vBp,S) -  logproposal(vSwp, vBp,  S_prop ) - logc;
        if alpARc > log(rand)
            flag = 1;
        end
        if or(count0>1000,notfound==1)
            break
        end
    end
    if flag == 1
        % MH-step
        alpAR = logtarget(vB,S) -  logproposal(vSw, vB,  S_prop ) - logc;
        if alpAR < 0
            alpMH = 0;
        elseif alpARc < 0
            alpMH = - alpAR;
        else
            alpMH = alpARc - alpAR;
        end
        if alpMH > log(rand)
            vSw = vSwp;
            vB = vBp;
            acc_rate =   acc_rate + 1;
        end
    end
    
    B = reshape(S_b'*vB,ntil,ntil);
    Sig = B*B'; iSig = inv(Sig);
    
    %% Save draws
    alpha_draws(:,irep ) = vec(Alpha);
    vecB_draws(:,irep ) = vec(B.*sign(diag(B))');
end
alpha_draws = alpha_draws(:,input.nburn+1:end);
vecB_draws = vecB_draws(:,input.nburn+1:end);
acc_rate = acc_rate/irep;

logtarget = @(xvecB) lpriorB(xvecB) ;
if input.compute_prior == 1
    %% Draw B | alpha, Y  by Accept Reject - Metropolis Hastings (AR-MH) step
    S_prop = S + S0aug;  % Update Proposal distribution for B
    if irep<input.nburn % Tune the constant during the burn-in...
        logc = (irep-1)/irep*logc + 1/irep*( logtarget(vB,S) - logproposal(vSw, vB,  S_prop )) ;
    end
    flag = 0;
    count0 = 0;
    while flag == 0
        count0 = count0 + 1;
        [vSwp, vBp, notfound ] = proposal_draw(S0aug);
        alpARc = logtarget(vBp,S) -  logproposal(vSwp, vBp,  S0aug ) - logc;
        if alpARc > log(rand)
            flag = 1;
        end
        if or(count0>1000,notfound==1)
            break
        end
    end
    if flag == 1
        % MH-step
        alpAR = logtarget(vB,S) -  logproposal(vSw, vB,  S_prop ) - logc;
        if alpAR < 0
            alpMH = 0;
        elseif alpARc < 0
            alpMH = - alpAR;
        else
            alpMH = alpARc - alpAR;
        end
        if alpMH > log(rand)
            vSw = vSwp;
            vB = vBp;
            acc_rate =   acc_rate + 1;
        end
    end
    
    B = reshape(S_b'*vB,ntil,ntil);
    
end











 
% compute the Bayes Factor: 
logdenominator = zeros(1e04,2);
for irep = 1:1e04
    Vphi_prior = kron(kappa_phi*eye(n),iwishrnd(S_eta,v_eta));  
    mu_prior = zeros(n*k,1);
    % Relevance:  
    S_r = eye(length(vec(Phi_r))); S_f = S_r;
    S_r(vec(Phi_r)==0,:) = []; 
    S_f(vec(Phi_r)==1,:) = [];
    mu_r = S_r*mu_prior + (S_r*Vphi_prior*S_f')*((S_f*Vphi_prior*S_f')\(S_f*vec(Phi_i) - S_f*mu_prior));
    V_r = S_r*Vphi_prior*S_r' - (S_r*Vphi_prior*S_f')*((S_f*Vphi_prior*S_f')\(S_r*Vphi_prior*S_f')');
    phi0_r = zeros(length(mu_r),1);
    logdenominator(irep,1) = - 0.5*length(mu_r)*log(2*pi) - 0.5*LogAbsDet(V_r)-0.5*( (phi0_r-mu_r)'*(V_r\(phi0_r-mu_r)) );  
    % Exogeneity:  
    S_r2 = eye(length(vec(Phi_r))); S_f2 = S_r2;
    S_r2(vec(Phi_r)==1,:) = []; 
    S_f2(vec(Phi_r)==0,:) = [];
    mu_r2 = S_r2*mu_prior + (S_r2*Vphi_prior*S_f2')*((S_f2*Vphi_prior*S_f2')\(S_f2*vec(Phi_i) - S_f2*mu_prior));
    V_r2 = S_r2*Vphi_prior*S_r2' - (S_r2*Vphi_prior*S_f2')*((S_f2*Vphi_prior*S_f2')\(S_r2*Vphi_prior*S_f2')');
    phi0_r2 = zeros(length(mu_r2),1); 
    logdenominator(irep,2) = - 0.5*length(mu_r2)*log(2*pi) - 0.5*LogAbsDet(V_r2)-0.5*( (phi0_r2-mu_r2)'*(V_r2\(phi0_r2-mu_r2)) );
end    

%logBFs = lognominator - logdenominator;
maxlogdenom = max(logdenominator);
logdenominatorhat = maxlogdenom+log(mean(exp(logdenominator-maxlogdenom)));
logBFs = (lognominator - logdenominatorhat);
maxlognom = max(lognominator);
lognominatorhat = maxlognom + log(mean(exp(lognominator-maxlognom)));
BFs = exp(lognominatorhat - logdenominatorhat);
scatter(lognominator(:,2),lognominator(:,1))

[maxnom,maxidx]=max(lognominator(:,2)); 
Bmax = reshape(vecB_draws(:,maxidx),ntil,ntil);
% 
% idx_s = logBFs(:,2)>0;
% Bset = reshape(vecB_draws(:,idx_s),ntil,ntil,sum(idx_s));
% for i = 1:size(Bset,3) 
%     Bset(:,input.targetshock,i) = Bset(:,input.targetshock,i).*sign(Bset(n,input.targetshock,i));
% end
% Bimp = squeeze(Bset(:,input.targetshock,:));
% 

output.BFs = BFs;
output.logBFs = logBFs; 
output.alphas = alpha_draws;
output.vecBs = vecB_draws;
output.input = input;
output.S_b = S_b;
output.Ika = Ika;
output.accs = acc_rate;

disp('check')
end


%% Some functions needed during the algorithm are put here:
function lprior = lnpriorB( vB, S_b, S0, v0 )
% Prior for B
ntil = size(S0,2);
Btil = reshape(S_b'*vB,ntil,ntil);
lprior = - (v0+ntil)*LogAbsDet(Btil) - .5*trace( S0*inv(Btil*Btil') );
end


function lprior = lnpriorBnew2( vB, S_b, S0, v0 , S_rphi)
% Prior for B
ntil = size(S0,2);
n = size(S_rphi,1); 
k = ntil - n;

S_11 = S0(1:n,1:n); S_12 = S0(1:n,n+1:n+k);
S_22 = S0(n+1:n+k,n+1:n+k);
S_22dot1 = S_22 - S_12'*(S_11\S_12);
Btil = reshape(S_b'*vB,ntil,ntil);
B = Btil(1:n,1:n);
Sig_eta = Btil(n+1:end,n+1:end)*Btil(n+1:end,n+1:end)';
vPhi = S_rphi*vec(Btil(n+1:end,1:n)); 
pm = S_rphi*vec(S_12'*(S_11\B));
pV = S_rphi*kron(B'*(S_11\B),  Sig_eta )*S_rphi';
%pV = S_rphi*kron( inv(S_11) , B'*Sig_eta*B)*S_rphi';
dPhim = vPhi-pm;
lprior1 = - (v0+n)/2*LogAbsDet(B*B') - .5*trace( S_11*inv(B*B') );
lprior2 = - (v0+n+k)/2*LogAbsDet(Sig_eta) - .5*trace( S_22dot1*inv(Sig_eta) );
lprior3 = -0.5*LogAbsDet(pV) - 0.5*(dPhim'*(pV\dPhim));
lprior = lprior1 + lprior2 + lprior3;
end



function lprior = lnpriorBnew( vB, S_b, S0, v0 , kappa_phi, S_eta, v_eta ,S_rphi)
% Prior for B
n = size(S0,2);
k = size(S_eta,2); 
ntil = n + k;
Btil = reshape(S_b'*vB,ntil,ntil);
B = Btil(1:n,1:n); 
vPhi = S_rphi*vec(Btil(n+1:end,1:n)); 
Sig_eta = Btil(n+1:end,n+1:end)*Btil(n+1:end,n+1:end)';
Sig_phi = S_rphi*kron(kappa_phi*eye(n),Sig_eta)*S_rphi';  
lprior = - (v0+n)*LogAbsDet(B) - .5*trace( S0*inv(B*B') )...
    - (v_eta+k)/2*LogAbsDet(Sig_eta) - .5*trace( S_eta*inv(Sig_eta) )...
    -0.5*LogAbsDet(Sig_phi) - 0.5*(vPhi'*(Sig_phi\vPhi));
end

function llh = llh_B(vB,S,S_b,T,ntil)
% Likelihood for B
B = reshape(S_b'*vB,ntil,ntil);
llh = - T*LogAbsDet(B) -0.5*trace((B\S)/B');
end
function logpdf = pstar_eval(vSw, vBp, Dnn, S, T, f, ZA, S_b)
% Density proposal draw
ntil = size(S,2);
Sigma = reshape(Dnn*vSw(1:size(Dnn,2)),ntil,ntil);
Dfx = Gradp(f,vBp);
if sum(sum(isnan(ZA)))==ntil^2
    log_volume = 0.5*LogAbsDet(Dfx'*Dfx);
else
    fz = @(xvB)zero_restr(xvB,ZA,S_b);
    Dhx = Gradp(fz,vBp);
    N = Dfx*null(Dhx);
    log_volume = 0.5*LogAbsDet(N'*N);
end
% log_volume = 0;
logwish = logiwishpdf(Sigma, S ,T);
logpdf = logwish + log_volume; % only up to normalizing constant of iWish pdf !
end


function zr = zero_restr(vB,ZA,S_b)
ntil = size(ZA,2);
B = reshape(S_b'*vB,ntil,ntil);
A = inv(B);
zr = A(ZA==1);
end



function lpdf = logiwishpdf(Sigma, S, v)
d = size(Sigma,2);
lpdf =  -0.5*(v+d+1)*LogAbsDet(Sigma)-0.5*trace(S/Sigma);
end

function [vSw,vBp, notfound,count] = pstar_draw(S, T, ntil, k, ztil, stil, Z1, Z2, S1, S2, Dnn, S_b, W, RInstr)
% Generates a draw from the AR-MH proposal s.t. sign and zero restrictions
n = ntil-k;
flag = 0;
count = 0;
notfound = 0;
while flag==0
    Sig = iwishrnd(S,T); %inv(wish(inv(S),T));
    w = [Draw_w(n,ztil(1:n)); Draw_w(k,ztil(n+1:end))];
    vSw = [ vech(Sig); w];
    vBp = f_map(vSw, ntil, k, ztil,Z1, Z2, Dnn, S_b, W);
    % vSw2 = f_map_inv(vBp, ntil, k, ztil, Z1, Z2, S_b, W);
    % vBp2 = f_map(vSw2, ntil, k, ztil,Z1, Z2, Dnn, S_b, W);
    B = reshape(S_b'*vBp,ntil,ntil);
    flag = 1;
    if nargin>13
        if RInstr(B)==1
        else
            flag = 0;
            continue
        end
    end
    G = [inv(B),B']';
    for j = find(stil~=0)'
        if j<=n
            if sum(S1{j}*G(:,j)>0)==length(S1{j}*G(:,j))
            elseif sum(S1{j}*G(:,j)>0)==0
            else
                flag = 0;
                continue
            end
        else
            if sum(S2{j-k}*G(:,j)>0)==length(S2{j-k}*G(:,j))
            elseif sum(S2{j-k}*G(:,j)>0)==0
            else
                flag = 0;
                continue
            end
        end
    end
    count = count + 1;
    if count>4*1e5
        notfound = 1;
        flag = 1;
    end
end


end


function w = Draw_w(n,ztil)
w = zeros(sum(1:n)-sum(ztil),1);
idx = 0 ;
for j = 1:n
    s = n + 1 - j - ztil(j);
    x1j = randn(s,1);
    w(idx+1:idx+s) = (x1j./norm(x1j));
    idx = idx + s;
end
end




function vecB = f_map(x, ntil, k, ztil, Z1, Z2, Dnn, S_b, W)
% maps [vech(Sigma), w] to S_B*vec(B)]
n = ntil-k;
nS = size(Dnn,2);
Sigma = reshape(Dnn*x(1:nS),ntil,ntil);
P = chol(Sigma)';
G = [inv(P),P']';
w = x(nS+1:end);
Q1 = zeros(n,n);
idx = 0;
for j = 1:n
    s = n + 1 - j - ztil(j) ;
    wj = w(idx+1:idx+s);
    Mj = [Q1(:,1:j-1)'; Z1{j}*G(:,1:n); W{j,1}];
    [K,R]=qr(Mj');
    for i = n-s+1:n
        if R(i,i)<0
            K(:,i)= -K(:,i);
        end
    end
    Kj = K(:,n-s+1:n);
    Q1(:,j) = Kj*wj;
    idx = idx + s;
end
Q2 = zeros(k,k);
for j = 1:k
    s = k + 1 - j - ztil(n+j);
    Mj = [Q2(:,1:j-1)'; Z2{j}*G(:,n+1:end); W{j,2}];
    wj = w(idx+1:idx+s);
    [K,R]=qr(Mj');
    for i = k-s+1:k
        if R(i,i)<0
            K(:,i)=-K(:,i);
        end
    end
    Kj = K(:,k-s+1:k);
    Q2(:,j) = Kj*wj;
    idx = idx + s;
end
Q = eye(ntil);
Q(1:n,1:n) = Q1;
Q(n+1:end,n+1:end) = Q2;
B = P*Q;
vecB = S_b*vec(B);
end

function vPcholw = f_map_inv(vB, ntil, k, ztil, Z1, Z2, S_b, W)
%vecB = f_g(x, ntil, k, ztil, Z1, Z2, Dnn, S_b, W)
%F_H Summary of this function goes here
%   Detailed explanation goes here
n = ntil-k;
B = reshape(S_b'*vB,ntil,ntil);
[Qp,R]=qr(B');
for i = 1:size(R,2)
    if R(i,i)<0
        Qp(:,i) = - Qp(:,i);
    end
end
Pchol = B*Qp;
Q1 = Qp(1:n,1:n)';
Q2 = Qp(n+1:end,n+1:end)';
G = [inv(Pchol),Pchol']';
w = [];
idx = 0;
for j = 1:n
    s = n + 1 - j - ztil(j) ;
    Mj = [Q1(:,1:j-1)'; Z1{j}*G(:,1:n); W{j,1}];
    [K,R]=qr(Mj');
    for i = n-s+1:n
        if R(i,i)<0
            K(:,i)= -K(:,i);
        end
    end
    Kj = K(:,n-s+1:n);
    wj = Kj\Q1(:,j);
    w = [w;wj];
    idx = idx + s;
end
for j = 1:k
    s = k + 1 - j - ztil(n+j);
    Mj = [Q2(:,1:j-1)'; Z2{j}*G(:,n+1:end); W{j,2}];
    [K,R]=qr(Mj');
    for i = k-s+1:k
        if R(i,i)<0
            K(:,i)=-K(:,i);
        end
    end
    Kj = K(:,k-s+1:k);
    wj = Kj\Q2(:,j);
    w = [w; wj];
    idx = idx + s;
end

vPcholw = [ vech(Pchol*Pchol');  w];
end

