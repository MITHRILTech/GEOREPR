function [cHF, cH, cF] = eqHF1(T, unit, cHF0, cH0, cF0, I)

T(T<25) = 25; I(I<0) = 0;
Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;
rho = -0.0025*TdegC.^2-0.1249*TdegC+1005.2;
e = -1.602 * 10^-19; % C, Electron Charge
NA = 6.022 * 10^23; % Avogadro’s number
kB = 1.38 * 10^-23; %m2kgs-2K-1, Boltzmann constant
Ew = relativepermittivityWater(T, unit); % wiki
Ew(Ew<0) = 0;
E = 8.8542E-12; % Fm-1, permittivity of free space
A = e^2 * (2*e^2*NA*rho./(E*Ew*kB.*Tab)).^0.5 ./ (2.302585*8*pi*E*Ew*kB.*Tab); %
B = (2*e^2*NA*rho./(E*Ew*kB.*Tab)).^0.5 * 10^-9; %
A(A==inf) = 0; B(B==inf) = 0;

gammaH = 10 .^ (-A * 1^2 .* I.^0.5 ./ (1 + B * 9.0 .* I.^0.5));
gammaHF = ones(size(T));
gammaF = 10 .^ (-A * 2^2 .* I.^0.5 ./ (1 + B * 0.35 .* I.^0.5));
gammaHF2 = gammaF;

K1 = 10 .^ (-2.75 - 0.014*Tab - 295./Tab + 1.91*log10(Tab));
K2 = 10 .^ (-4.389 + 0.00229*Tab + 77.62./Tab);

aHF0 = cHF0 .* gammaHF; aa = aHF0;
aH0 = cH0 .* gammaH; bb = aH0;
aF0 = cF0 .* gammaF; cc = aF0;

K0 = bb .* cc ./ aa;
K0(isnan(K0)) = 0; % FUCK 12.09am 02/04/2019
case1 = (K0 <= K1); % forward reaction
case2 = (K0 > K1); % backward reaction


a1 = 1;
b1 = (bb + cc + K1);
c1 = bb.*cc-K1.*aa;
x01 = (-b1 + (b1.^2 - 4*a1.*c1).^0.5) ./ (2 * a1);
x02 = (-b1 - (b1.^2 - 4*a1.*c1).^0.5) ./ (2 * a1);
signx11 = sign(x01) + sign(x02);
signCase02 = (signx11 == 2*ones(size(signx11)));
signCase00 = (signx11 == zeros(size(signx11)));
signCaseN01 = (signx11 == -1*ones(size(signx11)));
signCase01 = (signx11 == ones(size(signx11)));
x11 = (signCase02+signCase01) .* min(x01, x02) + (signCase00+signCaseN01) .* max(x01, x02);

a2 = 1;
b2 = -(bb + cc - K1);
c3 = bb.*cc-K1.*aa;
x1 = (-b2 + (b2.^2 - 4*a2.*c3).^0.5) ./ (2 * a2);
x2 = (-b2 - (b2.^2 - 4*a2.*c3).^0.5) ./ (2 * a2);
signx22 = sign(x1) + sign(x2);
signCase2 = (signx22 == 2*ones(size(signx22)));
signCase0 = (signx22 == zeros(size(signx22)));
signCaseN1 = (signx22 == -1*ones(size(signx22)));
signCase1 = (signx22 == ones(size(signx22)));
x22 = (signCase2+signCase1) .* min(x1, x2) + (signCase0+signCaseN1) .* max(x1, x2);

x = case1 .* x11 + case2 .* x22;

aHF = case1 .* (aa - x) + case2 .* (aa + x);
aH = case1 .* (bb + x) + case2 .* (bb - x);
aF = case1 .* (cc + x) + case2 .* (cc - x);

cHF = aHF ./ gammaHF;
cF = aF ./ gammaF;
cH = aH ./ gammaH;

end

