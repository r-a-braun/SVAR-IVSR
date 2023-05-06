 function output = SVAR_BB2021( input )
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

%% Set zero Constraints on vec(A) if exo==1
if input.exo == 1
    Alre = zeros(ntil);
    Alre(1:n,1:n) = ones(n);
    RestrLag = [ones(input.c,ntil);repmat(Alre,p,1)];
else
    RestrLag = [ones(input.c,ntil);repmat(ones(ntil),p,1)];
end
[na1,na2]=size(RestrLag);
Ika = speye(na1*na2);
Ika(vec(RestrLag)==0,:)=[];
nA = size(Ika,1);


%% Set Regressor Matrix 
lmX = lagmatrix(Ytilraw,1:p);
Xtil = [ones(T,input.c), lmX(p+1:end,:)];
Ytil = Ytilraw(p+1:end,:); 
if isfield(input,'n_T')==1 % Training sample exists:
    Y0 = Ytil(1:input.n_T,:);
    X0 = Xtil(1:input.n_T,:);
    u0 = zeros(input.n_T,n+k);
    
    for i = 1:n 
        if input.c==1
            if input.p>0
                Xi = X0(:,[1,1+i]); % just AR(1)
            else
                Xi = X0(:,[1]); % just AR(1)
            end
        else
            if input.p>0
                Xi = X0(:,i); % just AR(1)
            else
                Xi = []; % just AR(1)
            end 
        end
        bi = (Xi'*Xi)\(Xi'*Y0(:,i));
        u0(:,i) = Y0(:,i) - Xi*bi;
    end
    u0(:,n+1:end) = Y0(:,n+1:end);
    if input.priorind==1
       S0 = u0'*u0  ;  
       S0 = blkdiag(S0(1:n,1:n),S0(n+1:end,n+1:end));
    else
       S0  = u0'*u0;   
    end 
    v0 = input.n_T;   
    sig2vars = diag(S0/v0);
    kappa = [input.lam, input.lam/3, 1e04] ; 
    [mu_p,Var_p] = prior_Minnesota(n+k,p,input.prior_pers,input.c,kappa,sig2vars);
    a_0 = Ika*vec(mu_p);
    V_0 = sparse(Ika*diag(vec(Var_p))*Ika');
    iV_0 = inv(V_0); 
    iV0a0 = iV_0*a_0; 
    Xtil = Xtil(input.n_T+1:end,:);
    Ytil = Ytil(input.n_T+1:end,:); %Yexport = Ytil(input.n_T-input.p+1:end,:)
    Tp = Tp - input.n_T;
    T = Tp - p;
    
    
    %% Compute Bayes Factor as pretest of invertibility  
    nreg =   input.c;
    Omega_pr = zeros(ntil*p + nreg,1);
    Omega_pr(1:nreg) = kappa(3); % intercept
    for i=1:p
        Omega_pr((i-1)*ntil+nreg+1:nreg+i*ntil) =  (kappa(1))*(1/(i^2))./sig2vars;
    end
    % plot( kron(diag(diag(S0/v0)),diag(Omega))-diag(Var_p2(:)) ) % 
    kreg = size(Omega_pr,1); 
    iy1 = 1:n; % variables number 1 to 13 are in y1 
    iy2 = 1:ntil; iy2(iy1)=[];
    iy2X  = input.c + kron((0:(p-1))*ntil,ones(size(iy2))) + repmat(iy2,1,p);  
    Bprior = mu_p;  Sprior = S0; vprior = v0; Qprior =  sparse(1:kreg,1:kreg,Omega_pr);%diag(Omega); 
    logpriorpdf = matricvstudentpdf(zeros(length(iy2X),length(iy1)), ...
        Bprior(iy2X,iy1), Qprior(iy2X,iy2X)\eye(length(iy2X)), Sprior(iy1,iy1), vprior - length(iy2)); 
    % posterior pdf 
    iOmpost =  sparse(1:kreg,1:kreg,1./Omega_pr) + Xtil'*Xtil ; 
    Bst = iOmpost\(sparse(1:kreg,1:kreg,Omega_pr)\Bprior + Xtil'*Ytil); 
    Spost = Sprior + Bprior'*sparse(1:kreg,1:kreg,1./Omega_pr)*Bprior + Ytil'*Ytil  - Bst'*iOmpost*Bst; 
    vpost = vprior + length(Xtil);
    Qpost = inv(iOmpost);
    logpostpdf = matricvstudentpdf(zeros(length(iy2X),length(iy1)), ...
        Bst(iy2X,iy1), Qpost(iy2X,iy2X)\eye(length(iy2X)), Spost(iy1,iy1), vpost - length(iy2));  
    output.InvertibilityPreTest = 2*(logpriorpdf - logpostpdf);
     
    
    
    
    
    
    
    
else 
    S0 = eye(ntil); 
    v0 = ntil+2;  
    S0 = diag(v0*var_is);
    a_0 = zeros(nA,1); 
    V_0 = 100*eye(nA);
    iV_0 = inv(V_0); 
    iV0a0 = iV_0*a_0; 
end 
vYtil = vec(Ytil);
output.S0 = S0;
output.v0 = v0; 

%% Set identifying constraints on vec(B) and other functions
B_r = ones(ntil);
B_r(1:n,n+1:end) = 0;% thats the zero block we usually have 
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
%%  
lpriorB = @(vB)lnpriorB( vB, S_b, S0, v0 )  + input.pr_B(vB,S_b);
%lpriorB = @(vB)lnpriorBnew( vB, S_b, S0, v0, kappa_phi, S_eta, v_eta,S_rphi); % + input.logpvB2(vB)
% lpriorB = @(vB)lnpriorBnew2( vB, S_b, S0, v0 , S_rphi); % + input.logpvB2(vB)

%% Define some functions needed during the algorithm
v_prop = T + v0;
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
SSE = U'*U;
S_prop =  U'*U + S0;
ninit = 1000;
IS_w = zeros(1,ninit);
vBpcands = zeros(size(S_b,1),ninit);
for j = 1:ninit % Choose a high density point of this density
    [vSwp, vBp] = proposal_draw(S_prop);
    IS_w(j) =   logtarget(vBp,SSE) -  logproposal(vSwp,vBp,S_prop) ;  
    vBpcands(:,j) =  vBp;
end 
% var(IS_w(1:j))
[~,idx_max] =  max(IS_w);
vB = vBpcands(:,idx_max);
vSw = finv(vB);
B = reshape(S_b'*vB,ntil,ntil); 
Sig = B*B';

%% Save the following quantities:
alpha_draws = zeros(na1*na2,input.nrep);
vecB_draws = zeros(ntil^2,input.nrep);

%% MCMC Algorithm:
logc  = 0;
acc_rate = 1;
it_print = 1000; 
for irep = 1:input.nrep + input.nburn
    if mod(irep,it_print) == 0
        clc;
        disp(['Fraction of MCMC completed:',' ', num2str(irep/(input.nrep + input.nburn))])
        disp(['Acceptance rate:',' ', num2str(acc_rate/irep)]) 
    end 
    % Draw Alpha | Sigma, Y
    iSig = inv(B*B');
    iV_p = iV_0 + Ika*kron(iSig,Xtil'*Xtil)*Ika';
    p_m = iV_p\(iV0a0 + Ika*vec(Xtil'*Ytil*iSig));
    stat = 1; countstat = 1;
    while and(stat==1,countstat<1000)
        alphac = p_m + chol(iV_p,'lower')'\randn(nA,1);
        Alphac = reshape(Ika'*alphac,na1,ntil);
        if p > 0 
            stat = StationarityCheck2(Alphac,ntil,p);
            if stat == 0
                Alpha = Alphac; alpha = alphac;
            end
        else
            stat = 0;
        end
        countstat = countstat + 1;
    end
    S = (Ytil - Xtil*Alpha)'*(Ytil - Xtil*Alpha);
    
    %% Draw B | alpha, Y  by Accept Reject - Metropolis Hastings (AR-MH) step
    S_prop = S + S0;  % Update Proposal distribution for B
    if irep<input.nburn % Tune the constant during the burn-in...
        logc = (irep-1)/irep*logc + 1/irep*( logtarget(vB,S) - logproposal(vSw, vB,  S_prop )  + log(2)) ;
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
    if irep>input.nburn 
        alpha_draws(:,irep-input.nburn ) = vec(Alpha);
        vecB_draws(:,irep-input.nburn ) = vec(B);
    end 
end
 

acc_rate = acc_rate/irep;

if input.draw_prior == 1 
    % Generate also draws from the prior of B (v_0, S_0)
    % compute the Bayes Factor:
    vecB_draws_prior = zeros(ntil^2,input.nrep);
    acc_rate_pr = 0;
    if isfield(input,'SRinstru')==1 % Training sample exists:
        proposal_draw = @(xS)pstar_draw(S0, v0, ntil, k, ztil, stil, input.Z1, input.Z2, input.S1, input.S2, Dnn, S_b, W, RInstr);
    else
        proposal_draw = @(xS)pstar_draw(S0, v0, ntil, k, ztil, stil, input.Z1, input.Z2, input.S1, input.S2, Dnn, S_b, W);
    end
    [vSw, vB] = proposal_draw(S0);
    logproposal = @(xvPw ,xvB )pstar_eval(xvPw, xvB, Dnn, S0, v0, finv, input.ZR(:,:,1), S_b) ;
    logtarget = @(xvecB ) lpriorB(xvecB) ;  % not normalized 
    logc = 0;
    for irep = 1:input.nrep+input.nburn % Choose a high density point of this density
        if mod(irep,it_print) == 0
            clc;
            disp([irep/(input.nrep + input.nburn),acc_rate_pr/irep]);
        end
        
        if irep<input.nburn % Tune the constant during the burn-in...
            logc = (irep-1)/irep*logc + 1/irep*( logtarget(vB) - logproposal(vSw, vB)  + log(2)) ;
        end
        flag = 0;
        count0 = 0;
        while flag == 0
            count0 = count0 + 1;
            [vSwp, vBp, notfound ] = proposal_draw(S0);
            alpARc = logtarget(vBp) -  logproposal(vSwp, vBp) - logc;
            if alpARc > log(rand)
                flag = 1;
            end
            if or(count0>1000,notfound==1)
                break
            end
        end
        if flag == 1
            % MH-step
            alpAR = logtarget(vB) -  logproposal(vSw, vB) - logc;
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
                acc_rate_pr =   acc_rate_pr + 1;
            end
        end 
        
        if irep>input.nburn 
            vecB_draws_prior(:,irep-input.nburn ) = vec(reshape(S_b'*vB,ntil,ntil));
        end
    end 
    output.vecBs_pr = vecB_draws_prior;
    output.accs_pr = acc_rate_pr/(input.nrep+input.nburn);
end
 
output.alphas = alpha_draws;
output.vecBs = vecB_draws; 
output.input = input;
output.S_b = S_b;
output.Ika = Ika;
output.accs = acc_rate;

%disp('check')
end


%% Some functions needed during the algorithm are put here:
function lprior = lnpriorB( vB, S_b, S0, v0 )
% Prior for B
ntil = size(S0,2);
Btil = reshape(S_b'*vB,ntil,ntil);
lprior = - (v0+ntil)*LogAbsDet(Btil) - .5*trace( S0*inv(Btil*Btil') );
end 

% 
% function lprior = lnpriorBnew( vB, S_b, S0, v0 , kappa_phi, S_eta, v_eta ,S_rphi)
% % Prior for B
% n = size(S0,2);
% k = size(S_eta,2);
% ntil = n + k;
% Btil = reshape(S_b'*vB,ntil,ntil);
% B = Btil(1:n,1:n);
% vPhi = S_rphi*vec(Btil(n+1:end,1:n));
% Sig_eta = Btil(n+1:end,n+1:end)*Btil(n+1:end,n+1:end)';
% Sig_phi = S_rphi*kron(kappa_phi*eye(n),Sig_eta)*S_rphi';
% lprior = - (v0+n)*LogAbsDet(B) - .5*trace( S0*inv(B*B') )...
%     - (v_eta+k)/2*LogAbsDet(Sig_eta) - .5*trace( S_eta*inv(Sig_eta) )...
%     -0.5*LogAbsDet(Sig_phi) - 0.5*(vPhi'*(Sig_phi\vPhi));
% end








