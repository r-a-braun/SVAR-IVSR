function [ Z1,S1, Z2,S2 ] = get_Z_S( ZR,SR,n,k )
%GET_Z_S Summary of this function goes here
%   Detailed explanation goes here
%% Identification restrictions:
ntil = n + k; 
%%% Construct the matrices...
NSs = 2 ; % number of objects in F(A_{0},A_{+}) to which we impose sign restrictios: F(THETA)=[A_{0};L_{0},,...}] (
NSz = 2; %number of objects in F(A_{0},A_{+}) to which we impose zero restrictios: F(THETA)=[A_{0};L_{0};L_{inf}]
for i=1:n+k
    if i<=n
        ZRi = squeeze(ZR(:,i,:));
        idxi = find(isnan(vec(ZRi))==0);
        Z1{i} = zeros(length(idxi),ntil*NSz);
        for ii = 1:length(idxi)
            Z1{i}(ii,idxi(ii)) = 1;
        end
        SRi =  squeeze(SR.SIGN(:,i,:)) ; 
        vSRi = vec(SRi);
        idxi = find(isnan(vSRi)==0)';
        S1{i} = zeros(length(idxi), ntil*NSs);
        for ii = 1:length(idxi)
            S1{i}(ii,idxi(ii)) = vSRi( idxi(ii) );
        end
    else
        ZRi = squeeze(ZR(:,i,:));
        idxi = find(isnan(vec(ZRi))==0);
        Z2{i-n} = zeros(length(idxi),ntil*NSz);
        for ii = 1:length(idxi)
            Z2{i-n}(ii,idxi(ii)) = 1;
        end
        SRi =  squeeze(SR.SIGN(:,i,:)) ; 
        vSRi = vec(SRi);
        idxi = find(isnan(vSRi)==0)';
        S2{i-n} = zeros(length(idxi),ntil*NSs);
        for ii = 1:length(idxi)
            S2{i-n}(ii,idxi(ii)) = vSRi(idxi(ii));
        end
    end
end

end

