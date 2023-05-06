

function logpdf = pstar_eval(vSw, vBp, Dnn, S, T, f, ZA, S_b)
% Density proposal draw
ntil = size(S,2);
Sigma = reshape(Dnn*vSw(1:size(Dnn,2)),ntil,ntil);
Dfx = Gradp(f,vBp);
%Dfx = jacobianest(f, vBp);

if sum(sum(isnan(ZA)))==ntil^2
    log_volume = 0.5*LogAbsDet(Dfx'*Dfx); 
else
    fz = @(xvB)zero_restr(xvB,ZA,S_b);
    Dhx = Gradp(fz,vBp);
%     B = reshape(S_b'*vBp,ntil,ntil);
%     A = inv(B);
%     zr = A(ZA==1);
%     -kron(B',B)
   % Dhx = NumericalDerivative(fz, vBp);
    %Dhx = jacobianest(fz, vBp);
    N = Dfx*null(Dhx);
    log_volume = 0.5*LogAbsDet(N'*N);
end
% log_volume = 0;
logwish = logiwishpdf(Sigma, S ,T);
logpdf = logwish + log_volume; % only up to normalizing constant of iWish pdf !

% 
% 
% f_h_inv = @(vBp) [vec(reshape(S_b'*vBp,ntil,ntil)*reshape(S_b'*vBp,ntil,ntil)') ;...
%     vec(chol(reshape(S_b'*vBp,ntil,ntil)*reshape(S_b'*vBp,ntil,ntil)','lower')\reshape(S_b'*vBp,ntil,ntil))]; 
% %vSvQ = f_h_inv(A(:));
% % disp(Sig-reshape(vSvQ(1:n^2),n,n))  % check mapping
% % disp(Q-reshape(vSvQ(n^2+1:end),n,n))  % check mapping
%  
% % Evaluate Log Volume Element numerical and as in proposition 1
% Df = Gradp(f_h_inv, vBp); 
% lve_numerical = 0.5*LogAbsDet(Df'*Df);  
end