function [SDS, Astotal] = sinkDS(cDSinj, cDSi, Astotali, T, unit, pH, cNaCl, dt)

cDSi(isnan(cDSi)) = 0;
pH(isnan(pH)) = 7;
cNaCl(isnan(cNaCl)) = 0;

% clear all;
% [cDSinj, cDSi, Astotali, T, unit, pH, cNaCl, dt] = deal(1.04, 1.04, 0, 160, 'degC', 7, 0.09, 4671.1);

% [cDSinj, cDSi, Astotali, T, unit, pH, cNaCl, dt] = deal(1.0425, 0.6595, 2809, 160, 'degC', 7, 0.09, 36000);

%=============================================================%
%
% This function is used to model silica polymerisation in water by following Weres et al. (1980)
% at 50-100 degreeC (may be safely extrapolated at least up to 150 degreeC) and pH<8
%
% cSAppmi: initial dissolved silica concentration in ppm
% TdegC: temperature in degC
% pH: equivalent pH at room temperature
% cNaCl: ionic strength
% tmax: time
% indt: induction time switch
%
%=============================================================%

% Discretise time
Nt = 60; % number of nodes
Dt = dt / Nt;
t = 0 : Dt : dt;

% Astotal = zeros(length(t),1);
% SSI = zeros(length(t),1);
% n = zeros(length(t),1);
% Z = zeros(length(t),1);
% Acrit = zeros(length(t),1);
% dF = zeros(length(t),1);
% ff = zeros(length(t),1);
% Rmd = zeros(length(t),1);
% IN = zeros(length(t),1);
% N = zeros(length(t),1);
% As = zeros(length(t),1);
% cDS = zeros((length(t)+1),1);

% Constants
QLP = 3.34 * 10^25; % kg-1, Lothe-Pound factor
kBergs = 1.38 * 10^-16; % ergsK-1, Boltzmann constant
rhoNSA = 2.21 * 10^28; % m-3, The number density of SiO2 units in solids AS
rhoSilica = 2196; % kg/m3, Density of amorphous silica
h = 0.45;


Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;

fNaCl =10.^(-0.51*cNaCl.^0.5./(1+cNaCl.^0.5));
aNaCl = cNaCl.*fNaCl; % convert to activity

pHnom = pH + log10(aNaCl./0.069);

%Surface tension estimation
ipH = 10.^(-0.75924 + 0.58993*(pH-7.6) -0.11292 * (pH-7.6).^2);
ipHnom = 10.^(-0.75924 + 0.58993*(pHnom-7.6) -0.11292 * (pHnom-7.6).^2); %
I = 0.45 * ipH + 0.55 * ipHnom;
gammaergs = (63.68 - (0.049 + 0.2174 * I) .* Tab) * 10^4; % ergs m-2

% ASSUME: IDENTICAL AND FIXED
rho = -0.0025*TdegC.^2-0.1249*TdegC+1005.2; % kgm-3, density of water

St = 10.^(0.0977+75.84./Tab);

% f, f', p137
if pH<=5.97
    f0 = 10.^(pH-7.6);
    fpH = 10.^(f0./(1+6.2*f0));
    fpHprime = fpH/0.118913;
else 
    fpH = 10.^(pH - 7.6 - 2.113 * log10(1+10.^((pH-7.6)/2.113)) - (pH-7.6)./(9.6538 +...
        1.7901 * (pH-7.6) + 4.1811 * (pH-7.6).^2));
    fpHprime = fpH/0.118913;
end

if pHnom<=5.97
    f0nom = 10.^(pHnom-7.6);
    fpHnom = 10.^(f0nom./(1+6.2*f0nom));
    fpHnomprime = fpHnom/0.118913;
else 
    fpHnom = 10.^(pHnom - 7.6 - 2.113 * log10(1+10.^((pHnom-7.6)/2.113)) - (pHnom-...
        7.6)./(9.6538 + 1.7901 * (pHnom-7.6) + 4.1811 * (pHnom-7.6).^2));
    fpHnomprime = fpHnom/0.118913;
end

F = h * fpHprime + (1 - h) * fpHnomprime;
kOH = 10.^(3.1171 - 4296.6./Tab);

IS = cNaCl;
Ce = solubilityAS(T, unit, pH, IS);

SSI0 = cDSinj./Ce;
SSIi = cDSi./Ce;

% if SSI0 > ones(size(SSI0)) && SSIi > ones(size(SSI0))
%     
%     rcrit = 2 * gammaergs ./ (rhoNSA * kBergs .* Tab .* log(SSI0));
% else
%     rcrit = 0;
% end

testSSI1 = (SSI0 > ones(size(SSI0)));
testSSI2 = (SSIi > ones(size(SSI0)));
testSSI = (0.5 * (testSSI1 + testSSI2) >= ones(size(testSSI1)));
rcrit = testSSI .*...
    (2 * gammaergs ./ (rhoNSA * kBergs .* Tab .* log(SSI0)));




for i = 1 : length(t)
    
    cDS(:,:,1) = cDSi;
    SSI(:,:,i) = cDS(:,:,i)./Ce;

    % Homogeneous nucleation, Weres et al 1980
    n(:,:,i) = 4 * pi * rhoNSA * rcrit.^3 / 3;
    Z(:,:,i) = (2/3) * (3./(4*pi*rhoNSA*n(:,:,i).^2)).^(1/3) .* (gammaergs./(kBergs * Tab)).^(1/2);
    Acrit(:,:,i) = 4 * pi * rcrit.^2;
    dF(:,:,i) = gammaergs .* Acrit(:,:,i) / 3 ./ (kBergs * Tab);
    
    % Find the molecular deposition rate
    if SSI(:,:,i) > St
        ff(:,:,i) = St.^5 + 5 * St.^4 * (SSI(:,:,i) - St);
    else
        ff(:,:,i) = SSI(:,:,i).^5;
    end
    Rmd(:,:,i) = F .* kOH .* ff(:,:,i) .* (1-1./SSI(:,:,i)) *0.001/0.01/0.01/60;
    % kg SiO2 m-2 s-1, The molecular deposition rate
    
    % Find the nucleation rate
    IN(:,:,i) = QLP * Z(:,:,i) .* Rmd(:,:,i) / rhoSilica * rhoNSA .* Acrit(:,:,i) .* ...
        exp(-dF(:,:,i)) .* rho; % particles / m3 / s
    N(:,:,i) = IN(:,:,i) * Dt; % particles / m3
    As(:,:,i) = N(:,:,i) * 4 * pi .* rcrit.^2; % m2 / m3
    
    if i==1
        Astotal(:,:,i) = Astotali + As(:,:,i);
    else
        Astotal(:,:,i)  = Astotal(:,:,i-1) + As(:,:,i);
    end
    
    if ((Rmd(:,:,i) .* Astotal(:,:,i) * Dt) < (cDS(:,:,i) - Ce))
        cDS(:,:,i+1) = cDS(:,:,i) - Rmd(:,:,i) * Dt .* Astotal(:,:,i);
    elseif Rmd(:,:,i) > 0
        cDS(:,:,i+1) = Ce;
    else
        cDS(:,:,i+1) = cDS(:,:,i);
    end
    
    testM1 = (((Rmd(:,:,i) .* Astotal(:,:,i) * Dt) < (cDS(:,:,i) - Ce)));
    testM2 = (Rmd(:,:,i) > zeros(size(Rmd(:,:,i))));
    testM = testM1 + testM2;
    
    case1 = (testM == 2*ones(size(testM)));
    case2 = (testM == ones(size(testM)));
    case3 = (testM == zeros(size(testM)));
    
    cDS(:,:,i+1) = case1 .* (cDS(:,:,i) - Rmd(:,:,i) * Dt .* Astotal(:,:,i)) + case2 .* Ce + case3 .* cDS(:,:,i);
    
end

if dt == 0
    SDS = zeros(size(T));
    Astotal = zeros(size(T));
else
    SDS = (cDS(:,:,length(t)) - cDSi) / dt;
    SDS(isnan(SDS)) = 0;
    Astotal = Astotal(:,:,length(t));
    Astotal(isnan(Astotal)) = 0;
end

end


