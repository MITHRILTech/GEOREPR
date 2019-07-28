close all; clear all; clc;
% SILNUC is the work of Weres et al. (1980) 
% which models the kinetics of silica polymerisation.
% The semi-emprical formulas were fitted from observations at temperatures of 50-100 degreeC 
% can be safely extrapolated at least up to 150 degreeC, 
% The model should not be relied on very much when pH>8
% Default sample: temperature = 50degC, initial concentration of dissolved silica = 1000ppm, 
% initial pH = 6.23, initial concentration of NaCl = 0.09M.
% The following code is the reproduced MATLAB version of SILNUC.
% The possible coding errors or incorrect reproducing of SILNUC are entirely mine.
% NOTE: Estimated induction time is OVERRIDDEN. To enable it, one can
% comment out Line 118

%% Inputs
pH0 = 6.23; TdegC0 = 50; cDSppmi0 = 1000; cNaCl0 = 0.09;
tmax0 = timeSI(100, 'min'); % Time scale of interest

% To conveniently input serveral sets of data:
i = 1;
[pH, cNaCl, TdegC, cDSppmi, tmax] = deal(pH0(i), cNaCl0(i), TdegC0(i), cDSppmi0(i), tmax0(i));

%% Constants
QLP = 3.34 * 10^25; % kg-1, Lothe-Pound factor
kB = 1.38 * 10 ^ -23; % m2kgs-2K-1, Boltzmann constant
kBcal = 3.2976230 * 10^-24; % calK-1, Boltzmann constant
kBergs = 1.38 * 10^-16; % ergsK-1, Boltzmann constant
rhoNDS = 2.21 * 10^28; %m-3, The number density of SiO2 units in solids AS
rhoSilica = 2196; %kg/m3, Density of amorphous silica
molarmassSilica = 0.06008; %kg mol-1
molarmassSA = 96.113 * 10^-3; %kg mol-1
h = 0.45; B = 5.2; pK = 6.4;
NA = 6.022*10^23; % mol-1, Avogadro?s number
Rg = 8.3144598*10^7; %ergs K-1 mol-1, gas constant
SiO = 1.74 * 10^-7; %mol/m2, [SiO-] SURFACE DENSITY of dissociated silanol groups at pH=7.31 (Fleming 1985)
vSIO2 = 0.00925*10^-3; %m3/mol, molar volume of silicon dioxide, ref: https://www.azom.com/properties.aspx?ArticleID=1114
vSA = 54.3*10^-6; %m3/mol, molar volume of silicic acid, ref: http://www.chemspider.com/Chemical-Structure.14236.html

%% INPUT 1: Temperature
Tab = TdegC + 273; % convert to absolute temperature

%% INPUT 2: Background eletrolyte, i.e. concentration of NaCl
% cNaCl = 0.07; % molL-1
fNaCl =10^(-0.51*cNaCl^0.5/(1+cNaCl^0.5));
aNaCl = cNaCl*fNaCl; % convert to activity

%% INPUT 3: pH
% Equivalent pH, which treats the effects of H+ and NaCl similarly
pHnom = pH + log10(aNaCl/0.069);
% Surface tension estimation
ipH = 10^(-0.75924 + 0.58993*(pH-7.6) -0.11292 * (pH-7.6)^2);
ipHnom = 10^(-0.75924 + 0.58993*(pHnom-7.6) -0.11292 * (pHnom-7.6)^2); %
I = 0.45 * ipH + 0.55 * ipHnom;
gammaergs = (63.68 - (0.049 + 0.2174 * I) * Tab) * 10^4; % ergs m-2
gamma = gammaergs*0.2389*0.0000001; %cal m-2

%% INPUT 4: Initial monomeric silica concentration
% Nil

%% INPUT 5: Time discretisation
Nt = 60; % Number of timesteps: try increasing this if the outputs look incorrect.
dt =tmax/Nt*ones(Nt,1); t = 0:dt:tmax; tplot = 0:dt:(tmax-dt);

%% Pre-processing
rho = -0.0025*TdegC^2-0.1249*TdegC+1005.2; % kgm-3
mu = 0.000111; %Pas, Dynamic viscosity of injected and native reservoir fluids, ASSUME: IDENTICAL AND FIXED
nu = mu / rho; %m2s-1, Kinematic viscosity of injected and native reservoir fluids, ASSUME: IDENTICAL AND FIXED
rSA = (vSA / NA * (3/4) / pi) ^ (1/3); %averaged radius of Si(OH)4
DSA = kB * Tab / (3 * pi * 2 * rSA * mu); %m2s-1, Diffusion coefficient of silicic acid molecules

%% Matrix Initialisation
SSI = zeros(length(t),1); n = zeros(length(t),1); IN = zeros(length(t),1);
Z = zeros(length(t),1); Acrit = zeros(length(t),1); dF = zeros(length(t),1);
ff = zeros(length(t),1); Rmd = zeros(length(t),1); RmdSI = zeros(length(t),1);
cDS = zeros(length(t),1); cDSppm = zeros(length(t),1); cPS = zeros(length(t),1);
cPSppm = zeros(length(t),1); N = zeros(length(t),1); Ntotal = zeros(length(t),1);
As = zeros(length(t),1); Astotal = zeros(length(t),1); ta = zeros(length(t),1);

%% Solubility of dissolved silica
cTppm = 10 ^ (-731 / Tab+4.52);
K1 = 10 ^ (-2549 / Tab - 15.36 * 10^-6 * Tab^2);
betaCe = 10 ^ (- 0.51 * cNaCl^0.5 / (1 + 10^-10 * 0.32 * 10^8 * cNaCl^0.5));
Ceppm = cTppm * (1 + (10^pH * K1 / betaCe)); % Solubility in ppm
Ce = Ceppm / 1000 / 1000 * rho; % Solubility in kg/m3
SSIi = cDSppmi/1000/1000*rho/Ce; % Initial SSI
rcrit = 2 * gammaergs / (rhoNDS * kBergs * Tab * log(SSIi)); % Critical nucleation size
ni = 4 * pi * rhoNDS * rcrit^3 / 3;
Zi = (2/3) * (3/(4*pi*rhoNDS*ni^2))^(1/3) * (gammaergs/(kBergs * Tab))^(1/2);
Acriti = 4 * pi * rcrit^2;
St = 10^(0.0977+75.84/Tab);

% Efffect of pH and pHnom f, f', p137
if pH<=5.97
    f0 = 10^(pH-7.6); fpH = (f0/(1+6.2*f0)); fpHprime = fpH/0.118913;
else
    fpH = 10^(pH - 7.6 - 2.113 * log10(1+10^((pH-7.6)/2.113)) - (pH-7.6)/(9.6538 +...
        1.7901 * (pH-7.6) + 4.1811 * (pH-7.6)^2));
    fpHprime = fpH/0.118913;
end
if pHnom<=5.97
    f0nom = 10^(pHnom-7.6); fpHnom = (f0nom/(1+6.2*f0nom));
    fpHnomprime = fpHnom/0.118913;
else
    fpHnom = 10^(pHnom - 7.6 - 2.113 * log10(1+10^((pHnom-7.6)/2.113)) - (pHnom-...
        7.6)/(9.6538 + 1.7901 * (pHnom-7.6) + 4.1811 * (pHnom-7.6)^2));
    fpHnomprime = fpHnom/0.118913;
end
F = h * fpHprime + (1 - h) * fpHnomprime;
kOH = 10^(3.1171 - 4296.6/Tab);
dFi = gammaergs * Acriti / 3 / (kBergs * Tab);
if SSIi > St
    ffi = St^5 + 5 * St^4 * (SSIi - St);
else
    ffi = SSIi^5;
end
Rmdori = F * kOH * ffi * (1-1/SSIi);
% Estimation of induction time
tind = 1.08 * 10^-6 / Rmdori * (Zi * QLP * exp(-dFi) * (rcrit*100)^2 ) ^ (-0.25) *60; %s
tind = 0; %Induction time is OVERRIDEN here

for i = 1 : (length(t)-1)
    cDSppm(1) = cDSppmi; %ppm, Concentration of silicic acid in injected fluids
    cDS(1) = cDSppm(1)/1000/1000*rho; %kg/m3, ppm to kg/m3
    SSI(i) = cDS(i) / Ce;
    % Homogeneous nucleation
    rcrit(i) = 2 * gammaergs / (rhoNDS * kBergs * Tab * log(SSI(1)));
    n(i) = 4 * pi * rhoNDS * rcrit(i)^3 / 3;
    Z(i) = (2/3) * (3/(4*pi*rhoNDS*n(i)^2))^(1/3) * (gammaergs/(kBergs * Tab))^(1/2);
    Acrit(i) = 4 * pi * rcrit(i)^2;
    dF(i) = gammaergs * Acrit(i) / 3 / (kBergs * Tab);
    % Find the molecular deposition rate
    if SSI(i) > St
        ff(i) = St^5 + 5 * St^4 * (SSI(i) - St);
    else
        ff(i) = SSI(i)^5;
    end
    % The molecular deposition rate, in kg SiO2 m-2 s-1.
    RmdSI(i) = F * kOH * ff(i) * (1-1/SSI(i)) * 0.001/0.01/0.01/60;
    % The molecular deposition rate, in g SiO2 cm-2 min-1.
    Rmd(i) = F * kOH * ff(i) * (1-1/SSI(i));
    
    if t(i) < tind % Effect of induction time is taken into account here
        cDS(i+1) = cDS(i);
    else % after induction:
        % Find the nucleation rate
        IN(i) = QLP * Z(i) * RmdSI(i) / rhoSilica  * rhoNDS * Acrit(i) * ...
            exp(-dF(i))  * rho; % particles / m3 / s
        N(i) = IN(i) * dt(i); % Formed nuclei in the present timestep
        if i==1 % Total particle number
            Ntotal(i) = N(i);
        else
            Ntotal(i) = Ntotal(i-1) + N(i);
        end
        % Surface area of the formed nuclei in the present timestep
        As(i) = N(i) * 4 * pi * rcrit(i)^2; 
        % Total surface area of the formed nuclei 
        if i==1
            Astotal(i) = As(i);
        else
            Astotal(i)  = Astotal(i-1) + As(i);
        end
        cDS(i+1) = cDS(i) - RmdSI(i) * dt(i) * Astotal(i);
    end
    
    % Concentrations of dissoved(cDS) and polymerised(PS) sillica in the next timestep
    cDSppm(i+1) = cDS(i+1)*1000*1000/rho; % in ppm
    cPS(i+1) = cDS(1) - cDS(i+1); % in kg/m3
    cPSppm(i+1) = cPS(i+1)*1000*1000/rho; % in ppm
end

figure (1)
plot (t, cDSppm)
grid on
xlabel('Time (s)','fontsize',12,'fontweight','light');
ylabel('Concentration of dissolved silica (ppm)','fontsize',12,'fontweight','light');


