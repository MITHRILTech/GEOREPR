function rcrit = computeRcrit(cDS, T, unit, pH, IS)

Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;
kBergs = 1.38 * 10^-16; % ergsK-1, Boltzmann constant
rhoNSA = 2.21 * 10^28; % m-3, The number density of SiO2 units in solids AS
cNaCl = IS;

fNaCl =10.^(-0.51*cNaCl.^0.5./(1+cNaCl.^0.5));

aNaCl = cNaCl.*fNaCl; % convert to activity

pHnom = pH + log10(aNaCl/0.069);

ipH = 10.^(-0.75924 + 0.58993*(pH-7.6) -0.11292 * (pH-7.6).^2);
ipHnom = 10.^(-0.75924 + 0.58993*(pHnom-7.6) -0.11292 * (pHnom-7.6).^2); %

I = 0.45 * ipH + 0.55 * ipHnom;
gammaergs = (63.68 - (0.049 + 0.2174 * I) .* Tab) * 10^4; % ergs m-2

rho = -0.0025*TdegC.^2-0.1249*TdegC+1005.2; % kgm-3, density of water

% cDS = cDSppm / 1000 / 1000 * rho;

% Ceppm = solubilityAS(T, unit, pH, IS)*1000*1000./rho;
Ce = solubilityAS(T, unit, pH, IS);

SSI = cDS ./ Ce;

rcrit = 2 * gammaergs ./ (rhoNSA * kBergs * Tab .* log(SSI));
