function [alpha,beta]=GetBetaParameters(meanB,stdB)
% This function finds the parameter alpha and beta for the Beta
% distribution, such that the mode and the standard deviation of the
% distribution are equal to modeB and stdB, respectively
% Mode: (alpha-1)/(alpha+beta-2) for alpha and beta>1
% Mean: (alpha)/(alpha+beta) for alpha and beta>1
% Variance: (alpha*beta)/((alpha+beta)^2*(alpha+beta+1))
%
Bhat0=log(0.5)*ones(2,1);
% Optimisation:
options=optimset('MaxFunEvals',1000000);
Bhat=fminsearch('GetBetaParametersFM',Bhat0,options,meanB,stdB);
Bhat=exp(Bhat)+1;
alpha=Bhat(1);
beta=Bhat(2);