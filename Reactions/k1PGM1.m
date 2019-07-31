% SSI = [13.64 13.64 13.64 5.461 5.456 5.456 8.191 8.184 5.456];
% IS = [0.05 0.11 0.22 0.02 0.11 0.22 0.03 0.06 0.06];
% k1 = [0.0036 0.0039 0.0039 0.0027 0.0038 0.0044 0.002 0.0021 0.0015];

% SSI = [13.64 13.64 13.64 5.461 5.456 5.456];
% IS = [0.05 0.11 0.22 0.02 0.11 0.22];
% pHnom = [6.7669 7.0755 7.3407 6.399 7.0755 7.3407];
% k1 = [0.0036 0.0039 0.0039 0.0027 0.0038 0.0044];

% May only be safely applied when IS = 0.02 - 0.22 M, pH = 7, T = 25 degC

function k1 = k1PGM1(cDSinj, T, unit, pH, IS)



Tab = convert2SIunit(T, unit);
% TdegC = Tab - 273.15;

Ce = solubilityAS(T, unit, pH, IS); % kg m-3
SSI = cDSinj./Ce;
SSI(isnan(SSI)) = 0;

p00 =  -1.998e-05;
p10 =   0.0006216;
p01 =   0.0005915;
p20 =  -2.916e-05;
p11 =     0.00313;
p02 =   -0.001103;
p30 =  -1.162e-08;
p21 =  -0.0001616;
p12 =   -0.002896;
p03 =   -0.006444;

x = SSI;
y = IS;

if SSI <= 1
    k1 = 0;
else
    k1 = p00 + p10*x + p01*y + p20*x.^2 + p11*x.*y +...
        p02*y.^2 + p30*x.^3 + p21*x.^2.*y + p12*x.*y.^2 + p03*y.^3;
end

R = 1.9872036; %cal/k/mol, gas constant
E = 13.1 * 10^3; % cal/mol, activition energy

gamma25 = computeGamma(25, 'degC', pH, IS);
gamma = computeGamma(T, unit, pH, IS);
Ce25 = solubilityAS(25, 'degC', pH, IS);


k1 = gamma.*Ce./Tab * 298.15./(gamma25.*Ce25) .* k1 .* exp(- E/R * (1./Tab - 1/298.15));

end


% Linear model Poly12:
%      f(x,y) = p00 + p10*x + p01*y + p11*x*y + p02*y^2
% Coefficients (with 95% confidence bounds):
%        p00 =     0.00189  (0.001308, 0.002471)
%        p10 =     9.2e-05  (2.926e-05, 0.0001547)
%        p01 =     0.02012  (0.01209, 0.02815)
%        p11 =  -0.0007007  (-0.001126, -0.0002753)
%        p02 =    -0.03252  (-0.06416, -0.0008684)
% 
% Goodness of fit:
%   SSE: 5.792e-10
%   R-square: 0.9996
%   Adjusted R-square: 0.9982
%   RMSE: 2.407e-05

