function gamma = computeGamma(T, unit, pH, cNaCl)

Tab = convert2SIunit(T, unit);
% TdegC = Tab - 273.15;

fNaCl =10.^(-0.51*cNaCl.^0.5./(1+cNaCl.^0.5));
aNaCl = cNaCl.*fNaCl; % convert to activity

pHnom = pH + log10(aNaCl/0.069);

%Surface tension estimation
ipH = 10.^(-0.75924 + 0.58993*(pH-7.6) -0.11292 * (pH-7.6).^2);
ipHnom = 10.^(-0.75924 + 0.58993*(pHnom-7.6) -0.11292 * (pHnom-7.6).^2); %
I = 0.45 * ipH + 0.55 * ipHnom;
gammaergs = (63.68 - (0.049 + 0.2174 * I) .* Tab) * 10^4; % ergs m-2
gamma = gammaergs*0.2389*0.0000001; %cal m-2

end
