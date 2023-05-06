

function lpdf = logiwishpdf(Sigma, S, v)
d = size(Sigma,2);
lpdf =  -0.5*(v+d+1)*LogAbsDet(Sigma)-0.5*trace(S/Sigma);
end