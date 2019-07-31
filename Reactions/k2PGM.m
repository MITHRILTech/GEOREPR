% SSI = [13.64 13.64 13.64 5.461 5.456 5.456];
% IS = [0.05 0.11 0.22 0.02 0.11 0.22];
% pHnom = [6.7669 7.0755 7.3407 6.399 7.0755 7.3407];
% k1 = [0.0036 0.0039 0.0039 0.0027 0.0038 0.0044];
% k2 = [0.0107 0.0091 0.0056 0.005 0.0051 0.0065];

% Assume: k2 is a function of T, which has effects on D only, i.e. D2 =
% T2/T1*D1. Therefore, k2' = k2 * (T2 / T1)

function k2 = k2PGM(cDSinj, T, unit, pH, IS)

Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;

mu =  computeMu(T, unit); %Pas, Dynamic viscosity
mu25 =  computeMu(25, 'degC'); 

gamma25 = computeGamma(25, 'degC', pH, IS);
gamma = computeGamma(T, unit, pH, IS);

Ce = solubilityAS(T, unit, pH, IS); % kg m-3
Ce25 = solubilityAS(25, 'degC', pH, IS);
SSI = cDSinj ./ Ce;
SSI(isnan(SSI))=0;

p00 =  -0.0003503;
p10 =   0.0009594;
p01 =     0.02612;
p11 =   -0.004765;
p02 =     0.03119;

x = SSI;
y = IS;

if SSI <= 1
    k2 = 0;
else
    k2 = p00 + p10*x + p01*y + p11*x.*y + p02*y.^2;
end

k2 = gamma.*Ce./mu * mu25./(gamma25.*Ce25) .* k2 .* (Tab / 298.15);

end

% Linear model Poly12:
%      f(x,y) = p00 + p10*x + p01*y + p11*x*y + p02*y^2
% Coefficients (with 95% confidence bounds):
%        p00 =  -0.0003503  (-0.01188, 0.01118)
%        p10 =   0.0009594  (-0.0002837, 0.002203)
%        p01 =     0.02612  (-0.133, 0.1852)
%        p11 =   -0.004765  (-0.0132, 0.003664)
%        p02 =     0.03119  (-0.5959, 0.6583)

% 
% Goodness of fit:
%   SSE: 2.274e-07
%   R-square: 0.9919
%   Adjusted R-square: 0.9593
%   RMSE: 0.0004769