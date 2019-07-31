function Ceppm = solubilityASppm(T, unit, pH, IS)

T(T<25) = 25;
IS(IS<0) = 0;
% Comfirmed

Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;

rho = -0.0025*TdegC.^2-0.1249.*TdegC+1005.2;

cTppm = 10 .^ (-731 ./ Tab+4.52);
K1 = 10 .^ (-2549 ./ Tab - 15.36 * 10^-6 .* Tab .^2);
betaCe = 10 .^ (- 0.51 * IS .^0.5 ./ (1 + 10^-10 * 0.32 * 10^8 * IS .^0.5));
Ceppm = cTppm .* (1 + (10.^pH .* K1 ./ betaCe));

end
