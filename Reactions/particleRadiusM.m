function r = particleRadiusM(cDSinj, T, unit, pH, IS, tinj)

T(T<25) = 25;
IS(IS<0) = 0;

% M for matrix
% cDSinj, kg m-3, a matrix
% T, degC or K, a matrix
% pH, a matrix
% IS, M, a matrix

[~, Nz] = size(T);
tinj = [0; tinj; tinj(end)] * ones(1, Nz);

dmax = finalD(T, unit);

k1 = k1PGM1(cDSinj, T, unit, pH, IS);
k2 = k2PGM(cDSinj, T, unit, pH, IS);

rcritnm = 10^9*computeRcrit(cDSinj, T, unit, pH, IS);
rcritnm(rcritnm<0) = 0;


t1 = ((0.8*dmax*0.5).^2-rcritnm.^2)./k1;

t1(t1<0) = 0;


case1 = (tinj <= t1);
case2 = (tinj > t1);

r = (k1.*tinj + rcritnm.^2).^0.5 .* case1 + (k2.*(tinj-t1) + (0.8*dmax*0.5).^3).^(1/3) .* case2;
r(isnan(r))=0;

testM1 = (r<(dmax/2));
testM2 = (r>(dmax/2));

r = r .* testM1 + (dmax/2) .* testM2;



end


