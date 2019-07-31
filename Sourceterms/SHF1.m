function [SHF, SH, SF, SSiF6] = SHF1(T, unit, cHF0, cH0, cF0, cSiF60, I, m, a0, dt0)

cHF0(cHF0<0) = 0; cF0(cF0<0) = 0; cSiF60(cSiF60<0) = 0;

Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;
rho = -0.0025*TdegC.^2-0.1249*TdegC+1005.2;
Ar = 2 * [0; 2*pi.*m.cellcenters.x.*m.cellsize.x(2:end-1); 0]; 
H = a0/2*ones(size(m.cellsize.y)); 
V = Ar * H'; 
tmax = dt0;
tmin = 1;
Nt = 6;
qt = (tmax/tmin)^(1/(Nt-1));
for ii = 1:Nt
    if ii == 1
        t(ii) = tmin;
        dt(ii) = tmin;
    else
        t(ii) = tmin * qt^(ii-1);
        dt(ii) = t(ii) - t(ii-1);
    end
end

% Constant
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
gammaSiF6 = 10 .^ (-0.509 * 2^2 * I.^0.5 ./ (1 + (3.28*0.5*I.^0.5))); % Assume: di = 0.5
Kw = 10 .^ -(6E-5*25*ones(size(T)).^2-0.0286*25*ones(size(T))+14.463);
% Check if eq state is reached
cHdebug = cH0;
[cHF0, cH0, cF0] = eqHF1(T, unit, cHF0, cH0, cF0, I);
R = 8.314; % kJK-1mod-1

for i = 1:length(t)
    %% Initialisation
    if i == 1
        cHF = cHF0;
        cH = cH0;
        cF = cF0;
        cSiF6 = cSiF60;
    end
    
    aHF = cHF .* gammaHF;
    aH = cH .* gammaH;
    aF = cF .* gammaF;
    aSiF6 = cSiF6 .* gammaSiF6;
    pH = -log10(aH);
    cFmax = cHF + cF;
    cHFmax = cFmax;
    cHmax = cHF + cH;
    casepH1 = (pH > 3.16*ones(size(T)));
    casepH2 = (pH <= 3.16*ones(size(T)));
    Keq = casepH1 * 10^26.26 + casepH2 * 10^7.28;
    K = casepH1 .* (aSiF6 ./ (aH.^4 .* aF.^6)) + casepH2 .* (aSiF6 .* aH.^2 ./ aHF.^6); K(isnan(K)) = 0;
    caseJ11 = (K >= Keq); caseJ12 = (K < zeros(size(T))); caseJ13 = (cFmax <= zeros(size(T)));
    caseJ1a = caseJ11 + caseJ12 + caseJ13; caseJ1a(caseJ1a>1) = 1;
    caseJ1b = ones(size(T)) - caseJ1a;
    aHFlog = aHF; aHFlog(aHFlog<0) = 0; aHFlog(aHFlog==0) = 1; 
    
    J1 = 10 .^ (0.48 + 1.50 * log10(aHFlog) + 0.46 * pH - 1788.39 ./ Tab); J1 = caseJ1b .* J1; J1(isnan(J1)) = 0;
    
    J1(J1<0) = 0;
    mSiF6max = J1.* Ar .* dt(i); 
    mHmax = casepH1*4.*mSiF6max;
    mFmax = casepH1*6.*mSiF6max; 
    mHFmax = casepH2* 6.*mSiF6max;
   
    % logical matrix for pH > 3.16
    caseF1 = (cFmax.*V < mFmax); caseF1a =  ones(size(T)) - caseF1; % all available F 
    caseF1a = ones(size(T)) - caseF1;
    caseF2 = (cF.*V < mFmax); % all current F-
    caseF30 = caseF1 + caseF2; caseF30(caseF30>1) = 1;
    caseF3 = ones(size(T)) - caseF30; % Other cases
    caseH1 = (cHmax.*V < mHmax); caseH1a =  ones(size(T)) - caseH1; % all available F 
    caseH1a = ones(size(T)) - caseH1;
    caseH2 = (cH.*V < mHmax); % all current F-
    caseH30 = caseH1 + caseH2; caseH30(caseH30>1) = 1;
    caseH3 = ones(size(T)) - caseH30; % Other cases
    % logical matrix for pH < 3.16
    caseHF1 = (cFmax.*V < mHFmax); caseHF1a =  ones(size(T)) - caseHF1; % all available F 
    caseHF2 = (cHF.*V < mHFmax); % all current F-
    caseHF30 = caseHF1 + caseHF2; caseHF30(caseHF30>1) = 1;
    caseHF3 = ones(size(T)) - caseHF30; % Other cases
    % Product1
    cSiF6 = cSiF6 +...
        (caseF1 .* cFmax + caseF1a.*caseF2 .* cF + caseF3 .* ((mFmax./V))) / 6 +...
        (caseHF1 .* cHFmax + caseHF1a .* caseHF2 .* cHF +  caseHF3 .* (mHFmax./V)) / 6;
    cSiF6(isnan(cSiF6)) = 0;
    % Reactant1
    cF = caseF3 .* (cF - mFmax./V);
    cF(isnan(cF)) = 0;
    % Reactant1
    cH = caseH1 .* (Kw .^ 0.5) + caseH1a .* caseH2 .* (Kw .^ 0.5) +...
        caseH3 .* (cH - mHmax./V);
    cH(isnan(cH)) = 0;
    % Product2
    cH = cH + (caseHF3 .* (mHFmax./V)) / 3;
    cH(isnan(cH)) = 0;
    % Reactant2
    cHF = caseHF3 .* (cHF - mHFmax./V);
    cHF(isnan(cHF)) = 0;
    cF = cF - caseHF1 .* cF;
    
    [cHF, cH, cF] = eqHF1(T, unit, cHF, cH, cF, I);
    
    cFT = 6*cSiF6 + cHF + cF;
    
end

SHF =  (cHF - cHF0) / dt0;
SH =  (cH - cH0) / dt0;
SF =  (cF - cF0) / dt0;
SSiF6 =  (cSiF6 - cSiF60) / dt0;

end













