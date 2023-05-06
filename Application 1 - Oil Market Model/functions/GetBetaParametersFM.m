function FM=GetBetaParametersFM(x,meanB,stdB)
Bhat=exp(x)+1;
alpha=Bhat(1);
beta=Bhat(2);
varB=stdB^2;
FM=(meanB-(alpha)/(alpha+beta))^2+(varB-(alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))^2;




