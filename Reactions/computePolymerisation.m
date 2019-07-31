
function [cDSw, Astotal0] = computePolymerisation(cDSi, T, unit, pH, cNaCl, tmax, indt)

Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;

% clear all;
% [cDSi, TdegC, pH, cNaCl, tmax, indt] = deal(1.1,160,7,0.09,3600,0);

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
Nt = 10; % number of nodes
dt = tmax / Nt;
t = 0 : dt : tmax;

% Constants
QLP = 3.34 * 10^25; % kg-1, Lothe-Pound factor
kB = 1.38 * 10 ^ -23; % m2kgs-2K-1, Boltzmann constant
kBergs = 1.38 * 10^-16; % ergsK-1, Boltzmann constant
rhoNSA = 2.21 * 10^28; % m-3, The number density of SiO2 units in solids AS
rhoSilica = 2196; % kg/m3, Density of amorphous silica
h = 0.45;
NA = 6.022*10^23; % mol-1, Avogadro?s number
vSA = 54.3*10^-6; % m3/mol, molar volume of silicic acid,
% ref: http://www.chemspider.com/Chemical-Structure.14236.html

% Pre Processing
T = TdegC + 273; % convert to absolute temperature

fNaCl =10^(-0.51*cNaCl^0.5/(1+cNaCl^0.5));
aNaCl = cNaCl*fNaCl; % convert to activity
pHnom = pH + log10(aNaCl/0.069);

%Surface tension estimation
ipH = 10^(-0.75924 + 0.58993*(pH-7.6) -0.11292 * (pH-7.6)^2);
ipHnom = 10^(-0.75924 + 0.58993*(pHnom-7.6) -0.11292 * (pHnom-7.6)^2); %
I = 0.45 * ipH + 0.55 * ipHnom;
gammaergs = (63.68 - (0.049 + 0.2174 * I) * T) * 10^4;

% ASSUME: IDENTICAL AND FIXED
rho = -0.0025*TdegC^2-0.1249*TdegC+1005.2; % kgm-3, density of water

% Matrix Initialisation
SSI = zeros(length(t),1);
rcrit = zeros(length(t),1);
n = zeros(length(t),1);
Z = zeros(length(t),1);
Acrit = zeros(length(t),1);
dF = zeros(length(t),1);
ff = zeros(length(t),1);
Rmd = zeros(length(t),1);
Rmdor = zeros(length(t),1);
cDS = zeros(length(t),1);
cSAppm = zeros(length(t),1);
N = zeros(length(t),1);
Ntotal = zeros(length(t),1);
As = zeros(length(t),1);
Astotal = zeros(length(t),1);
IN = zeros(length(t),1);

St = 10^(0.0977+75.84/T);
% f, f', p137
if pH<=5.97
    f0 = 10^(pH-7.6);
    fpH = 10^(f0/(1+6.2*f0));
    fpHprime = fpH/0.118913;
else %if (pH0(m,i)<8.72915) && (pH0(m,i)>5.97)
    fpH = 10^(pH - 7.6 - 2.113 * log10(1+10^((pH-7.6)/2.113)) - (pH-7.6)/(9.6538 + 1.7901 * (pH-7.6) + 4.1811 * (pH-7.6)^2));
    fpHprime = fpH/0.118913;
end
if pHnom<=5.97
    f0nom = 10^(pHnom-7.6);
    fpHnom = 10^(f0nom/(1+6.2*f0nom));
    fpHnomprime = fpHnom/0.118913;
else %if (pH0(m,i)<8.72915) && (pH0(m,i)>5.97)
    fpHnom = 10^(pHnom - 7.6 - 2.113 * log10(1+10^((pHnom-7.6)/2.113)) - (pHnom-7.6)/(9.6538 + 1.7901 * (pHnom-7.6) + 4.1811 * (pHnom-7.6)^2));
    fpHnomprime = fpHnom/0.118913;
end
F = h * fpHprime + (1 - h) * fpHnomprime;
kOH = 10^(3.1171 - 4296.6/T);


%% SOLUBILITY
cTppm = 10 ^ (-731 / T+4.52);
K1 = 10 ^ (-2549 / T - 15.36 * 10^-6 * T^2);
betaCe = 10 ^ (- 0.51 * cNaCl^0.5 / (1 + 10^-10 * 0.32 * 10^8 * cNaCl^0.5));
Ceppm = cTppm * (1 + (10^pH * K1 / betaCe));
Ce = Ceppm / 1000 / 1000 * rho; %kg/m3

SSIi = cDSi/Ce;
rcriti = 2 * gammaergs / (rhoNSA * kBergs * T * log(SSIi));
ni = 4 * pi * rhoNSA * rcriti^3 / 3;
Zi = (2/3) * (3/(4*pi*rhoNSA*ni^2))^(1/3) * (gammaergs/(kBergs * T))^(1/2);
Acriti = 4 * pi * rcriti^2;
dFi = gammaergs * Acriti / 3 / (kBergs * T);
if SSIi > St
    ffi = St^5 + 5 * St^4 * (SSIi - St);
else
    ffi = SSIi^5;
end
Rmdori = F * kOH * ffi * (1-1/SSIi);
tind = 1.08 * 10^-6 / Rmdori * (Zi * QLP * exp(-dFi) * (rcriti*100)^2 ) ^ (-0.25) *60; %s

% Rmax = (3/4/(pi*rhoSilica) * (cSAppmi/1000/1000*rho - Ce))^(1/3);

% Induction time switch
if indt == 0
    tind = 0; %Induction time is OVERRIDEN
end


for i = 1 : (length(t)-1)
    
    cDS(1) = cDSi; %kg/m3, ppm to kg/m3
    SSI(i) = cDS(i) / Ce;
    
    % Homogeneous nucleation, Weres et al 1980
    rcrit(i) = 2 * gammaergs / (rhoNSA * kBergs * T * log(SSI(1)));
    n(i) = 4 * pi * rhoNSA * rcrit(i)^3 / 3;
    Z(i) = (2/3) * (3/(4*pi*rhoNSA*n(i)^2))^(1/3) * (gammaergs/(kBergs * T))^(1/2);
    Acrit(i) = 4 * pi * rcrit(i)^2;
    dF(i) = gammaergs * Acrit(i) / 3 / (kBergs * T);
    
    % Find the molecular deposition rate
    if SSI(i) > St
        ff(i) = St^5 + 5 * St^4 * (SSI(i) - St);
    else
        ff(i) = SSI(i)^5;
    end
    Rmd(i) = F * kOH * ff(i) * (1-1/SSI(i)) *0.001/0.01/0.01/60;
    % kg SiO2 m-2 s-1, The molecular deposition rate
    Rmdor(i) = F * kOH * ff(i) * (1-1/SSI(i));
    
    if t(i) < tind % Induction Time is taken into account here
        cDS(i+1) = cDS(i);
    else
        % Find the nucleation rate
        IN(i) = QLP * Z(i) * Rmd(i) / rhoSilica  * rhoNSA * Acrit(i) * ...
            exp(-dF(i)) * (1 - exp(-(t(i) + dt) * (4 * rhoNSA * Rmd(i) * Acrit(i)...
            * Z(i)^2))) * rho; % particles / m3 / s
        N(i) = IN(i) * dt;
        
        if i==1 % Total particle number
            Ntotal(i) = N(i);
        else
            Ntotal(i) = Ntotal(i-1) + N(i);
        end
        
        As(i) = N(i) * 4 * pi * rcrit(i)^2;
        
        if i==1
            Astotal(i) = As(i);
        else
            Astotal(i)  = Astotal(i-1) + As(i);
        end
        
        if ((Rmd(i) * Astotal(i) * dt) < (cDS(i) - Ce))
            cDS(i+1) = cDS(i) - Rmd(i) * dt * Astotal(i);
        elseif Rmd(i) > 0
            cDS(i+1) = Ce;
        else
            cDS(i+1) = cDS(i);
        end
        
    end
    
    %     cSAppm = cSA(i) * 1000 * 1000 / rho;
    cDSppm(i) = cDS(i) * 1000 * 1000 / rho;
    
end

cDSppm = [cDSppm cDSppm(i)];

% plot(t, cDSppm)

cDSw = cDS(length(t));
Astotal0 = Astotal(length(t)-1); 

end

