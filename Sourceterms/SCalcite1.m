function [SH, SH2SO4, SHSO4, SSO4, SCa, SCaHCO3, SH2CO3, SHCO3, SCO3] = SCalcite1(T, unit,...
    cNaCl0, cH0, cH2SO40, cHSO40, cSO40, cCa0, cCaHCO30, cCaT0, cH2CO30, cHCO30, cCO30, cCT0, m, a0, sCaCO3, dt0)

% clear all;
% close all;
% clc;
% [T, unit, cNaCl0, cH0, cH2SO40, cHSO40, cSO40, cCa0, cCaHCO30,...
%     cCaT0, cH2CO30, cHCO30, cCO30, cCT0, dx, dz, dt0] = deal(25, 'degC',...
%     0.09, 10^-7, 1.58e-5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1.5, 180);
cNaCl0(cNaCl0<0) = 0; cH0(cH0<0) = 0; cH2SO40(cH2SO40<0) = 0;
cHSO40(cHSO40<0) = 0; cSO40(cSO40<0) = 0; cCa0(cCa0<0) = 0;
cCaHCO30(cCaHCO30<0) = 0; cCaT0(cCaT0<0) = 0; cH2CO30(cH2CO30<0) = 0;
cHCO30(cHCO30<0) = 0; cCO30(cCO30<0) = 0; cCT0(cCT0<0) = 0;

Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;

rho = -0.0025*TdegC.^2-0.1249*TdegC+1005.2;
Ar = 2 * sCaCO3 .* [0; 2*pi.*m.cellcenters.x.*m.cellsize.x(2:end-1); 0]; % Calcite reaction surface area
% H = m.cellsize.y;
% H = m.cellcenters.y - m.cellcenters.y(Nzrock+1); H(1:Nzrock) = 0;
% H = [0; H; 0];%%%%%
% Hmax = max(m.cellcenters.y - m.cellcenters.y(Nzrock+1)); H = Hmax*ones(size(m.cellsize.y));
% H(1:Nzrock) = 0;
H = a0/2*ones(size(m.cellsize.y)); 
V = Ar * H'; % Assume: the calcite surface is a smooth flat plate

tmax = dt0;
tmin = 1;

% Nt = 60;
Nt = 6;
% dt = tmax / Nt; % ???
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

% Equilibrium constant
Kspcalcite = 10 .^ (10.22 - 0.0349*Tab - 2476./Tab); % activity product
Kspanhydrite = 10 .^ (4.287 - 0.02092*Tab - 423.7./Tab);
% Kspcalciteplot = Kspcalcite*ones(length(t)+1);
% Kspanhydriteplot = Kspanhydrite*ones(length(t)+1);
KaCaHCO3 = 10 .^ (1209.120 + 0.31294 * Tab - 34765.05./Tab - 478.782 * log10(Tab));
KdCaHCO3 = 1 ./ KaCaHCO3;
KdH2CO3 = 10 .^ (-356.3094 - 0.06091964*Tab + 21834.37./Tab + 126.8339*log10(Tab) - 1684915./Tab.^2);
KdHCO3 = 10 .^ (-107.8871 - 0.03252894 * Tab + 5151.79 ./ Tab + 38.92561*log10(Tab) - 563713.9 ./ Tab.^2);
KaCaCO3 = 10 .^ (-1228.732 - 0.299444 * Tab + 35512.75./Tab + 485.818 * log10(Tab));
KdCaCO3 = 1 ./ KaCaCO3;
% Kw = 10 .^ -(6E-5*TdegC.^2-0.0286*TdegC+14.463);
% Assume: logKw = -14, to find pH
Kw = 10 .^ -(6E-5*25*ones(size(T)).^2-0.0286*25*ones(size(T))+14.463);
KdHSO4 = K2H2SO4(Tab, unit);

% Rate constant
k1 = 10 .^ (0.198 - 444./Tab); % at stirring rate of about 1800 - 2300 rpm
k2 = 10 .^ (2.84 - 2177./Tab);

k3case1a = (TdegC <= 25*ones(size(TdegC)));
k3case1b = (TdegC >= 5*ones(size(TdegC)));
k3case1 = 0.5*(k3case1a + k3case1b);
k3case1(k3case1<1) = 0;
k3case2 = (TdegC > 25*ones(size(TdegC)));

k3 = k3case1 .* 10 .^ (-5.86 - 317./Tab) +...
    k3case2 .* 10 .^ (-1.10 - 1737./Tab);

% if TdegC <= 25 && TdegC >= 5
%     k3 = 10 ^ (-5.86 - 317/Tab);
% elseif TdegC > 25 %TdegC <= 48 && TdegC > 25
%     k3 = 10 ^ (-1.10 - 1737/Tab);
% end
%     k4 = KdHCO3 / KdCaCO3 * (k1 + 1/aH(i) * (k2 * aH2CO3(i) + k3 * aH2O));
k41 = k1 .* KdHCO3 ./ KdCaCO3;
k42 = k2 .* KdHCO3 ./ (KdCaCO3 .* KdH2CO3);
k43 = k3 .* KdHCO3 ./ (KdCaCO3 .* Kw);

% Eq. constants
K1 = Kspcalcite ./ KdHCO3;
K2 = Kspcalcite .* KdH2CO3 ./ KdHCO3;
K3 = Kspcalcite .* Kw ./ KdHCO3;

aH2O = ones(size(T)); % Assumption
% caseH2O0 = zeros(size(T));
caseH2O1 = zeros(size(T));

for i = 1:length(t)
    %% Initialisation
    if i == 1
        I(:,:,i) = cNaCl0;
        gammaH0 = 10 .^ (-A * 1^2 .* I(:,:,i).^0.5 ./ (1 + B * 9.0 .* I(:,:,i).^0.5));
        cHSO4(:,:,i) = cHSO40;
        cSO4(:,:,i) = cSO40;
        aH(:,:,i) = Kw .^ 0.5;
        cH(:,:,i) = aH(:,:,i) ./ 10 .^ (-A .* 1^2 .* I(:,:,i).^0.5 ./ (1 + B * 9.0 .* I(:,:,i).^0.5));
        aOH(:,:,i) = Kw ./ aH(:,:,i);
        cH2SO4(:,:,i) = cH2SO40;
        cST(:,:,i) = cH2SO4(:,:,i) + cHSO4(:,:,i) + cSO4(:,:,i);
        
        %         if cH2SO4(:,:,i) == 0
        [cHSO4(:,:,i+1), cSO4(:,:,i+1)] = eq1H2SO4(T, unit, cH(:,:,i), cST(:,:,i), I(:,:,i));
        %         else
        %             [cH2SO4(:,:,i), cH(:,:,i), cHSO4(:,:,i), cSO4(:,:,i)] =...
        %                 eqH2SO4(T, unit, cH(:,:,i), cH2SO4(:,:,i), cHSO4(:,:,i), cSO4(:,:,i), I(:,:,i));
        %         end
        
        
        aH(:,:,i) = cH(:,:,i) .* 10 .^ (-A * 1^2 .* I(:,:,i).^0.5 ./ (1 + B * 9.0 .* I(:,:,i).^0.5));
        
        caseI1 = (I(:,:,i) <= 0.1*ones(size(T)));
        caseI2 = (I(:,:,i) > 0.1*ones(size(T)));
        
        aSO4(:,:,i) = caseI1 .* (cSO4(:,:,i) .* 10 .^ (-0.509 * 1^2 * I(:,:,i).^0.5 ./ (1 + 3.28 * 0.4 * I(:,:,i).^0.5))) +...
            caseI2 .* (cSO4(:,:,i) .* 10 .^ (-0.509 * (I(:,:,i).^0.5 ./ (1 + I(:,:,i).^0.5) - 0.3 * I(:,:,i))));
        aHSO4(:,:,i) = caseI1 .* (cSO4(:,:,i) .* 10 .^ (-0.509 * (I(:,:,i).^0.5 ./ (1 + I(:,:,i).^0.5) - 0.3 * I(:,:,i)))) + ...
            caseI2 .* (cHSO4(:,:,i) .* 10 .^ (-0.509 * (I(:,:,i).^0.5 ./ (1 + I(:,:,i).^0.5) - 0.3 * I(:,:,i))));
        
        %         if I(:,:,i) <= 0.1
        %             aSO4(:,:,i) = cSO4(:,:,i) * 10 ^ (-0.509 * 1^2 * I(:,:,i)^0.5 / (1 + 3.28 * 0.4 * I(:,:,i)^0.5));
        %             aHSO4(:,:,i) = cHSO4(:,:,i) * 10 ^ (-0.509 * 1^2 * I(:,:,i)^0.5 / (1 + 3.28 * 0.4 * I(:,:,i)^0.5));
        %         else
        %             aSO4(:,:,i) = cSO4(:,:,i) * 10 ^ (-0.509 * (I(:,:,i)^0.5 / (1 + I(:,:,i)^0.5) - 0.3 * I(:,:,i)));
        %             aHSO4(:,:,i) = cHSO4(:,:,i) * 10 ^ (-0.509 * (I(:,:,i)^0.5 / (1 + I(:,:,i)^0.5) - 0.3 * I(:,:,i)));
        %         end
        
        cCT(:,:,i) = cCT0;
        cCaT(:,:,i) = cCaT0;
        
        [aH2CO3(:,:,i), aHCO3(:,:,i), aCO3(:,:,i), aCa(:,:,i), aCaCO30(:,:,i), aCaHCO3(:,:,i)] =...
            eq2CO2(T, unit, cCT(:,:,i), cCaT(:,:,i), aH(:,:,i));
        Kcalcite(:,:,i) = aCa(:,:,i).*aCO3(:,:,i);
        Kanhydrite(:,:,i) = aCa(:,:,i).*aSO4(:,:,i);
        pH(:,:,i) = -log10(aH(:,:,i));
        
    end
    
    %% Activity coeffcients
    gammaCa(:,:,i) = 10 .^ (-A * 2^2 .* I(:,:,i).^0.5 ./ (1 + B * 5.0 .* I(:,:,i).^0.5) + 0.165 * I(:,:,i));
    gammaHCO3(:,:,i) = 10 .^ (-A * 1^2 .* I(:,:,i).^0.5 ./ (1 + B * 5.4 .* I(:,:,i).^0.5) + 0.0 * I(:,:,i));
    gammaH2CO3(:,:,i) = ones(size(T));
    gammaCaHCO3(:,:,i) = gammaHCO3(:,:,i);
    gammaCO3(:,:,i) = 10 .^ (-0.509 * 2^2 * I(:,:,i).^0.5 ./ (1 + (3.28*0.5*I(:,:,i).^0.5))); % ???
    gammaCaCO3(:,:,i) = 10 .^ (-0.5 * I(:,:,i));
    gammaCaSO4(:,:,i) = 10 .^ (-0.45 * I(:,:,i));
    gammaH(:,:,i) = 10 .^ (-A * 1^2 .* I(:,:,i).^0.5 ./ (1 + B * 9.0 .* I(:,:,i).^0.5)); % 9 is the ion size parameter
    gammaOH(:,:,i) = 10 .^ (-A * 1^2 .* I(:,:,i).^0.5 ./ (1 + B * 3.5 .* I(:,:,i).^0.5)); % 3.5 is the ion size parameter
    gammaHSO4(:,:,i) = caseI1.*(10 .^ (-0.509 * 1^2 * I(:,:,i).^0.5 ./ (1 + 3.28 * 0.4 * I(:,:,i).^0.5))) +...
        caseI2.* (10 .^ (-0.509 * (I(:,:,i).^0.5 ./ (1 + I(:,:,i).^0.5) - 0.3 * I(:,:,i))));
    
    %     if I(i) <= 0.1
    %         gammaHSO4(i) = 10 ^ (-0.509 * 1^2 * I(i)^0.5 / (1 + 3.28 * 0.4 * I(i)^0.5));
    %     else
    %         gammaHSO4(i) = 10 ^ (-0.509 * (I(i)^0.5 / (1 + I(i)^0.5) - 0.3 * I(i)));
    %     end
    
    gammaSO4(:,:,i) = gammaHSO4(:,:,i);

    %% Rate constant
    aH(:,:,i) = cH(:,:,i) .* gammaH(:,:,i);
    cHmax(:,:,i) = cH(:,:,i) + cHSO4(:,:,i);
    aHmax(:,:,i) = cHmax(:,:,i) .* gammaH(:,:,i);
    
    caseJ11 = ((aCa(:,:,i) .* aHCO3(:,:,i) ./ aH(:,:,i) ./ K1) >= ones(size(T)));
    caseJ12 =  (aCa(:,:,i) .* aHCO3(:,:,i) ./ aH(:,:,i) ./ K1) < zeros(size(T));
    caseJ13 = (aHmax(:,:,i) < Kw.^0.5);
    caseJ14 = (Kcalcite(:,:,i) > Kspcalcite);
    caseJ1a = caseJ11 + caseJ12 + caseJ13+ caseJ14;
    caseJ1a(caseJ1a>1) = 1;
    caseJ1b = ones(size(T))  - caseJ1a;
    J1(:,:,i) = 0 .* caseJ1a + caseJ1b .*...
        (k1 .* aH(:,:,i) .* (1 - aCa(:,:,i) .* aHCO3(:,:,i) ./ aH(:,:,i) ./ K1));
    J1(isnan(J1)) = 0;
    %     if (1 - aCa(i) * aHCO3(i) / aH(i) / K1) <= 0 || (1 - aCa(i) * aHCO3(i) / aH(i) / K1) > 1 ||...
    %             aHmax(i) < Kw^0.5 || Kcalcite(i) > Kspcalcite
    %         J1(i) = 0;
    %     else
    %         J1(i) = k1 * aH(i) * (1 - aCa(i) * aHCO3(i) / aH(i) / K1); %molm-2s-1
    %     end
    
    caseJ21 = (aH2CO3(:,:,i) == zeros(size(T)));
    caseJ22 = (ones(size(T)) <= aCa(:,:,i) .* aHCO3(:,:,i) .^ 2 ./ aH2CO3(:,:,i) ./ K2);
    caseJ23 = ((aCa(:,:,i) .* aHCO3(:,:,i) .^ 2 ./ aH2CO3(:,:,i) ./ K2) < zeros(size(T)));
    caseJ24 = (Kcalcite(:,:,i) > Kspcalcite);
    caseJ2a = caseJ21 + caseJ22 + caseJ23+ caseJ24;
    caseJ2a(caseJ2a>1) = 1;
    caseJ2b = ones(size(T))  - caseJ2a;
    
    J2(:,:,i) = 0 .* caseJ2a + caseJ2b .*...
        (k2 .* aH2CO3(:,:,i) .* (1 - aCa(:,:,i) .* aHCO3(:,:,i) .^ 2 ./ aH2CO3(:,:,i) ./ K2));
    J2(isnan(J2)) = 0;
    
    %     if aH2CO3(i) == 0 || (1 - aCa(i) * aHCO3(i) ^ 2 / aH2CO3(i) / K2) <= 0 || (1 - aCa(i) *...
    %             aHCO3(i) ^ 2 / aH2CO3(i) / K2) > 1 || Kcalcite(i) > Kspcalcite
    %         J2(i) = 0;
    %     else
    %         J2(i) = k2 * aH2CO3(i) * (1 - aCa(i) * aHCO3(i) ^ 2 / aH2CO3(i) / K2); %molm-2s-1
    %     end
    
    caseJ31 = (aCa(:,:,i) .* aHCO3(:,:,i) .* aOH(:,:,i) ./ K3) >= zeros(size(T));
    caseJ32 = (aCa(:,:,i) .* aHCO3(:,:,i) .* aOH(:,:,i) ./ K3 < zeros(size(T)));
    caseJ33 = (Kcalcite(:,:,i) > Kspcalcite);
    caseJ3a = caseJ31 + caseJ32 + caseJ33;
    caseJ3a(caseJ3a>1) = 1;
    caseJ3b = ones(size(T))  - caseJ3a;
    J3(:,:,i) = 0 .* caseJ3a + caseJ3b .*...
        (k3 .* aH2O .* (1 - aCa(:,:,i) .* aHCO3(:,:,i) .* aOH(:,:,i) ./ K3));
%     J3(isnan(J3)) = 0;
    
    %     if (1 - aCa(i) * aHCO3(i) * aOH(i) / K3) <= 0 || (1 - aCa(i) * aHCO3(i) * aOH(i) / K3) > 1 ||...
    %             Kcalcite(i) > Kspcalcite
    %         J3(i) = 0;
    %     else
    %         J3(i) = k3 * aH2O * (1 - aCa(i) * aHCO3(i) * aOH(i) / K3); %molm-2s-1
    %     end
    
    %     J4(i) = -k4 * aCa(i) * aHCO3(i); %molm-2s-1
    J41(:,:,i) = -k41 .* aCa(:,:,i) .* aHCO3(:,:,i); %molm-2s-1
    J42(:,:,i) = -k42 .* aCa(:,:,i) .* aHCO3(:,:,i) .^ 2; %molm-2s-1
    J43(:,:,i) = -k43 .* aCa(:,:,i) .* aHCO3(:,:,i) .* aOH(:,:,i); %molm-2s-1
    
    %% H+ attack
    m1max(:,:,i) = J1(:,:,i) .* Ar .* dt(i);
    caseH1 = ((aHmax(:,:,i) .* V - m1max(:,:,i)) < Kw.^0.5.*V);
    caseH1a = ones(size(T)) - caseH1;
    caseH2 = ((aH(:,:,i) .* V - m1max(:,:,i)) < Kw.^0.5.*V);
    caseH30 = caseH1 + caseH2;
    caseH30(caseH30>1) = 1;
    caseH3 = ones(size(T)) - caseH30;
    
    KdHSO4prime = K2H2SO4prime(Tab, unit, I(:,:,i));
    
    cCaT(:,:,i+1) = caseH1 .* (cCaT(:,:,i) + aHmax(:,:,i))+...
        caseH1a .* caseH2 .* (cCaT(:,:,i) + m1max(:,:,i) ./ V)+...
        caseH3 .* (cCaT(:,:,i) + m1max(:,:,i) ./ V);
    cCaT(isnan(cCaT)) = 0;
    cCT(:,:,i+1) = caseH1 .* (cCT(:,:,i) + aHmax(:,:,i))+...
        caseH1a .* caseH2 .* (cCT(:,:,i) + m1max(:,:,i) ./ V)+...
        caseH3 .* (cCT(:,:,i) + m1max(:,:,i) ./ V);
    cCT(isnan(cCT)) = 0;
    
    aH(:,:,i+1) = caseH1 .* (Kw .^ 0.5)+...
        caseH1a .* caseH2 .* (Kw .^ 0.5)+...
        caseH3 .* ((aH(:,:,i) .* V -  m1max(:,:,i)) ./ V);
    
    aH(isnan(aH)) = 0;
    aHdebug = aH(:,:,i+1);
    
    [cHSO4(:,:,i+1), cSO4(:,:,i+1)] = eq1H2SO4(T, unit, cH(:,:,i), cST(:,:,i), I(:,:,i)); %?
    
    cSO4(:,:,i+1) = caseH1 .* (cHSO4(:,:,i) + cSO4(:,:,i))+...
        caseH1a .* caseH2 .* (cHSO4(:,:,i) - (m1max(:,:,i) - aH(:,:,i) .* V) + cSO4(:,:,i))+...
        caseH3 .* (cSO4(:,:,i+1));
%     cSO4(:,:,i+1) = caseH1 .* (cHSO4(:,:,i) + cSO4(:,:,i))+...
%         caseH1a .* caseH2 .* (cHSO4(:,:,i) - (m1max(:,:,i)./V - aH(:,:,i)) + cSO4(:,:,i))+...
%         caseH3 .* (cSO4(:,:,i+1));
    cHSO4(:,:,i+1) = 0*caseH1+...
        caseH1a .* caseH2 .* (cST(:,:,i) - cSO4(:,:,i))+...
        caseH3 .* (cHSO4(:,:,i+1));
    
    aSO4(:,:,i+1) = cSO4(:,:,i+1) .* gammaSO4(:,:,i);
    aHSO4(:,:,i+1) = cHSO4(:,:,i+1) .* gammaHSO4(:,:,i);
    
%     cH(:,:,i+1) = caseH1 .* (aH(:,:,i+1) ./ gammaH(:,:,i))+...
%         caseH1a .* caseH2 .* (KdHSO4prime .* cHSO4(:,:,i+1) ./ cSO4(:,:,i+1))+...
%         caseH3 .* (cH(:,:,i));
    cHa = caseH1 .* (aH(:,:,i+1) ./ gammaH(:,:,i));
    cHb = caseH1a .* caseH2 .* (KdHSO4prime .* cHSO4(:,:,i+1) ./ cSO4(:,:,i+1));
    cHb(isnan(cHb)) = 0; 
    cHc = caseH3 .* (cH(:,:,i));
    
    
    cH(:,:,i+1) = cHa + cHb + cHc;
    aH(:,:,i+1) = cH(:,:,i+1) .* gammaH(:,:,i);
    aHdebug = aH(:,:,i+1);
    cHdebug = cH(:,:,i+1);
    aOH(:,:,i+1) = Kw ./ aH(:,:,i+1);
    cOH(:,:,i+1) = aOH(:,:,i+1) ./ gammaOH(:,:,i);
    aOHdebug = aOH(:,:,i+1);
    cH2SO4(:,:,i+1) = zeros(size(T));
    
    
    %     if (aHmax(i) * V - m1max(i)) < Kw^0.5*V
    %         cCaT(i+1) = cCaT(i) + aHmax(i);
    %         cCT(i+1) = cCT(i) + aHmax(i);
    %         aH(i+1) = Kw ^ 0.5;
    %         cH(i+1) = aH(i+1) / gammaH(i);
    %         aOH(i+1) = Kw / aH(i+1);
    %         cSO4(i+1) = cHSO4(i) + cSO4(i);
    %         aSO4(i+1) = cSO4(i+1) * gammaSO4(i);
    %         cHSO4(i+1) = 0;
    %         aHSO4(i+1) = cHSO4(i+1) * gammaHSO4(i);
    %         cH2SO4(i+1) = 0;
    %     elseif (aH(i) * V - m1max(i)) < Kw^0.5*V
    %         cCaT(i+1) = cCaT(i) + m1max(i) / V;
    %         cCT(i+1) = cCT(i) + m1max(i) / V;
    %         aH(i) = Kw ^ 0.5;
    %         cSO4(i+1) = cHSO4(i) - (m1max(i) - aH(i) * V) + cSO4(i);
    %         aSO4(i+1) = cSO4(i+1) * gammaSO4(i);
    %         cHSO4(i+1) = cST - cSO4(i);
    %         aHSO4(i+1) = cHSO4(i+1) * gammaHSO4(i);
    %         KdHSO4prime = K2H2SO4prime(Tab, unit, I(i));
    %         cH(i+1) = KdHSO4prime * cHSO4(i+1) / cSO4(i+1);
    %         aH(i+1) = cH(i+1) * gammaH(i);
    %         aOH(i+1) = Kw / aH(i+1);
    %         cH2SO4(i+1) = 0;
    %     else
    %         cCaT(i+1) = cCaT(i) + m1max(i) / V;
    %         cCT(i+1) = cCT(i) + m1max(i) / V;
    %         aH(i) = (aH(i) * V -  m1max(i)) / V;
    %         cH(i) = aH(i) / gammaH(i);
    %         [cHSO4(i+1), cSO4(i+1)] = eq1H2SO4(T, unit, cH(i), cST(i), I(i));
    %         cH(i+1) = cH(i);
    %         cH2SO4(i+1) = cH2SO4(i);
    %         aH(i+1) = cH(i+1) * gammaH(i);
    %         aOH(i+1) = Kw / aH(i+1);
    %         aSO4(i+1) = cSO4(i+1) * gammaSO4(i);
    %         aHSO4(i+1) = cHSO4(i+1) * gammaHSO4(i);
    %     end
    
    %% H2CO3 attack
    m2max(:,:,i) = J2(:,:,i) .* Ar .* dt(i);
    caseH2CO301 = (aH2CO3(:,:,i) .* V < m2max(:,:,i));
    caseH2CO302 = ones(size(T)) - caseH2CO301;
    cCaT(:,:,i+1) = caseH2CO301.*(cCaT(:,:,i+1) + aH2CO3(:,:,i)) +...
        caseH2CO302.*(cCaT(:,:,i+1) + m2max(:,:,i) ./ V);
    cCaT(isnan(cCaT)) = 0;
    cCT(:,:,i+1) = caseH2CO301.*(cCT(:,:,i+1) + aH2CO3(:,:,i)) +...
        caseH2CO302.*(cCT(:,:,i+1) + m2max(:,:,i) ./ V);
    cCT(isnan(cCT)) = 0;
    %     if (aH2CO3(i) * V - m2max(i)) < 0
    %         cCaT(i+1) = cCaT(i+1) + aH2CO3(i);
    %         cCT(i+1) = cCT(i+1) + aH2CO3(i);
    %     else
    %         cCaT(i+1) = cCaT(i+1) + m2max(i) / V;
    %         cCT(i+1) = cCT(i+1) + m2max(i) / V;
    %     end
    
    %% H2O attack & Backward reaction
    
    
    
    for iii = 1:i
        
        caseH2O0 = (caseH2O1 == ones(size(T)));
        caseH2O1 = (Kcalcite(:,:,iii) > Kspcalcite);
        caseH2O1(caseH2O0==1) = 1;
        caseH2O2 = ones(size(T)) - caseH2O1;
        
        %         if Kcalcite(iii) > Kspcalcite
        %             m3max(i) = 0;
        %             m41max(i) = 0;
        %             m42max(i) = 0;
        %             m43max(i) = 0;
        %             break;
        %         else
        %             m3max(i) = J3(i) * Ar * dt(i);
        %             m41max(i) = J41(i) * Ar * dt(i);
        %             m42max(i) = J42(i) * Ar * dt(i);
        %             m43max(i) = J43(i) * Ar * dt(i);
        %         end
    end
    
    m3max(:,:,i) = 0.*caseH2O1 + caseH2O2.*(J3(:,:,i) .* Ar * dt(i));
    m41max(:,:,i) = 0.*caseH2O1 + caseH2O2.*(J41(:,:,i) .* Ar * dt(i));
    m42max(:,:,i) = 0.*caseH2O1 + caseH2O2.*(J42(:,:,i) .* Ar * dt(i));
    m43max(:,:,i) = 0.*caseH2O1 + caseH2O2.*(J43(:,:,i) .* Ar * dt(i));
    
    cCaT(:,:,i+1) = cCaT(:,:,i+1) + m3max(:,:,i) ./ V;
    cCaT(isnan(cCaT)) = 0;
    cCT(:,:,i+1) = cCT(:,:,i+1) + m3max(:,:,i) ./ V;
    cCT(isnan(cCT)) = 0;
    caseHOH1 = (aH(:,:,i+1) < Kw.^0.5);
    caseHOH2 = ((aH(:,:,i+1) - m3max(:,:,i) ./ V) <= Kw.^0.5);
    caseHOH30 = caseHOH1 + caseHOH2; caseHOH30(caseHOH30>1) = 1;
    caseHOH3 = ones(size(T)) - caseHOH30;
    
    aOH(:,:,i+1) = caseHOH1.* (aOH(:,:,i) + m3max(:,:,i) ./ V)+...
        caseHOH2.* (aOH(:,:,i) + (m3max(:,:,i) ./ V - aH(:,:,i+1) + Kw.^0.5)) +...
        caseHOH3.* (Kw ./ (aH(:,:,i+1) - m3max(:,:,i) ./ V));

    aOHdebug =  aOH(:,:,i+1);
    aOH(isnan(aOH)) = 0;
    
    aH(:,:,i+1) = caseHOH1.* (Kw ./ aOH(:,:,i+1))+...
        caseHOH2.* (Kw ./ aOH(:,:,i+1))+...
        caseHOH3.* (aH(:,:,i+1) - m3max(:,:,i) ./ V);
    
    aH(isnan(aH)) = 0;
    aHdebug = aH(:,:,i+1);
    
%     aOH(:,:,i+1) = caseHOH30.*aOH(:,:,i+1) + caseHOH3.*(Kw ./ aH(:,:,i+1));
%     aOH(isnan(aOH)) = 0;
    
    %     if aH(i+1) < Kw^0.5
    %         aOH(i+1) = aOH(i) + m3max(i) / V;
    %         aH(i+1) = Kw / aOH(i+1);
    %     elseif (aH(i+1) - m3max(i) / V) <= Kw^0.5
    %         aOH(i+1) = aOH(i) + (m3max(i) / V - aH(i+1) + Kw^0.5);
    %         aH(i+1) = Kw / aOH(i+1);
    %     else
    %         aH(i+1) = aH(i+1) - m3max(i) / V;
    %         aOH(i+1) = Kw / aH(i+1);
    %     end
    
    cH(:,:,i+1) = aH(:,:,i+1) ./ gammaH(:,:,i);
    cH(isnan(cH)) = 0;
    cHdebug = cH(:,:,i+1);
    cCaT(:,:,i+1) = cCaT(:,:,i+1) + m41max(:,:,i) ./ V + m42max(:,:,i) ./ V + m43max(:,:,i) ./ V;
    cCT(:,:,i+1) = cCT(:,:,i+1) + m41max(:,:,i) ./ V + m42max(:,:,i) ./ V + m43max(:,:,i) ./ V;
    cCaT(isnan(cCaT)) = 0; cCT(isnan(cCT)) = 0;
    
%     aH(:,:,i+1) = aH(:,:,i+1) - m41max(:,:,i) ./ V - m42max(:,:,i) ./ V; %???
    
    
    aOH(:,:,i+1) = Kw ./ aH(:,:,i+1); aH(isnan(aH)) = 0; aOH(isnan(aOH)) = 0;
    
    aHdebug = aH(:,:,i+1);
    aOHdebug =  aOH(:,:,i+1);
    [aH2CO3(:,:,i+1), aHCO3(:,:,i+1), aCO3(:,:,i+1), aCa(:,:,i+1), aCaCO30(:,:,i+1), aCaHCO3(:,:,i+1)] =...
        eq2CO2(T, unit, cCT(:,:,i+1), cCaT(:,:,i+1), aH(:,:,i+1));
    aH2CO3(isnan(aH2CO3)) = 0; aHCO3(isnan(aHCO3)) = 0; aCO3(isnan(aCO3)) = 0;
    aCa(isnan(aCa)) = 0; aCaCO30(isnan(aCaCO30)) = 0; aCaHCO3(isnan(aCaHCO3)) = 0;
    
    Kcalcite(:,:,i+1) = aCa(:,:,i+1).*aCO3(:,:,i+1);
    Kanhydrite(:,:,i+1) = aCa(:,:,i+1).*aSO4(:,:,i+1);
    
    %% Post
    % ??
    for j = 1 : 10
        
        %         jct(i) = j;
        
        casePost1 = (Kcalcite(:,:,i+1) > Kspcalcite);
        casePost2 = ((Kcalcite(:,:,i+1) - Kspcalcite) > 0.00001 .* Kspcalcite);
        casePost = casePost1 + casePost2;
        casePost(casePost > 1) = 1;
        casePost0 = ones(size(T)) - casePost;
        
        
        a = 1;
        b = -(aCa(:,:,i+1) + aCO3(:,:,i+1));
        c = aCa(:,:,i+1) .* aCO3(:,:,i+1) - Kspcalcite;
        x1 = (-b + (b.^2 - 4*a.*c).^0.5) ./ (2 * a);
        x2 = (-b - (b.^2 - 4*a.*c).^0.5) ./ (2 * a);
%         signx = sign(x1) .* sign(x2);
%         caseSamesign = (signx >= 0);
%         caseDiffsign = ones(size(signx)) - caseSamesign;
%         x = casePost .* (caseSamesign .* min(x1, x2) + caseDiffsign .* max(x1, x2)) +...
%             0 * casePost0;

        signx = sign(x1) + sign(x2);
        signCase2 = (signx == 2*ones(size(signx)));
        % signCaseN2 = (signx == -2*ones(size(signx)));
        signCase0 = (signx == zeros(size(signx)));
        signCaseN1 = (signx == -1*ones(size(signx)));
        signCase1 = (signx == ones(size(signx)));
        x = casePost .* ((signCase2+signCase1) .* min(x1, x2) + (signCase0+signCaseN1) .* max(x1, x2));
        
%         casex1 = (x1 >= zeros(size(T)));
%         casex2 = (x2 >= zeros(size(T)));
%         x1 = x1 .* casex1;
%         x2 = x2 .* casex2;
%         casex1min = (x1 <= x2);
%         casex2min = (x2 <= x1);
%         x = casePost .* (x1.*casex1min + x2.*casex2min) +...
%             0 * casePost0;
        
        aH2CO3(:,:,i+1) = x;
        
        cCT(:,:,i+1) = cCT(:,:,i+1) - x;
        cCaT(:,:,i+1) = cCaT(:,:,i+1) - x;
        
        [aH2CO3(:,:,i+1), aHCO3(:,:,i+1), aCO3(:,:,i+1), aCa(:,:,i+1), aCaCO30(:,:,i+1), aCaHCO3(:,:,i+1)] =...
            eq2CO2(T, unit, cCT(:,:,i+1), cCaT(:,:,i+1), aH(:,:,i+1));
        aH2CO3(isnan(aH2CO3)) = 0; aHCO3(isnan(aHCO3)) = 0; aCO3(isnan(aCO3)) = 0;
        aCa(isnan(aCa)) = 0; aCaCO30(isnan(aCaCO30)) = 0; aCaHCO3(isnan(aCaHCO3)) = 0;
        Kcalcite(:,:,i+1) = aCa(:,:,i+1).*aCO3(:,:,i+1);
        Kanhydrite(:,:,i+1) = aCa(:,:,i+1).*aSO4(:,:,i+1);
        
        %         if Kcalcite(i+1) > Kspcalcite && (Kcalcite(i+1) - Kspcalcite) > 0.00001 * Kspcalcite
        %
        %             a = 1;
        %             b = -(aCa(i+1) + aCO3(i+1));
        %             c = aCa(i+1) * aCO3(i+1) - Kspcalcite;
        %             x(1) = (-b + (b^2 - 4*a*c)^0.5) / (2 * a);
        %             x(2) = (-b - (b^2 - 4*a*c)^0.5) / (2 * a);
        %             x(x < 0) = [];
        %             aH2CO3 = min(x);
        %             cCT(i+1) = cCT(i+1) - min(x);
        %             cCaT(i+1) = cCaT(i+1) - min(x);
        %             [aH2CO3(i+1), aHCO3(i+1), aCO3(i+1), aCa(i+1), aCaCO30(i+1), aCaHCO3(i+1)] = eq2CO2(T, unit,...
        %                 cCT(i+1), cCaT(i+1), aH(i+1));
        %             Kcalcite(i+1) = aCa(i+1)*aCO3(i+1);
        %             Kanhydrite(i+1) = aCa(i+1)*aSO4(i+1);
        %
        %         else
        %             break
        %         end
    end
    
    %% Final
    
    cOH(:,:,i+1) = aOH(:,:,i+1) ./ gammaOH(:,:,i);
    cH2CO3(:,:,i+1) = aH2CO3(:,:,i+1) ./ gammaH2CO3(:,:,i);
    cHCO3(:,:,i+1) = aHCO3(:,:,i+1) ./ gammaHCO3(:,:,i);
    cCO3(:,:,i+1) = aCO3(:,:,i+1) ./ gammaCO3(:,:,i);
    cCa(:,:,i+1) = aCa(:,:,i+1) ./ gammaCa(:,:,i);
    cCaHCO3(:,:,i+1) = aCaHCO3(:,:,i+1) ./ gammaCaHCO3(:,:,i);
    cST(:,:,i+1) = cH2SO4(:,:,i+1) + cHSO4(:,:,i+1) + cSO4(:,:,i+1);
    
    
    cHdebug = cH(:,:,i+1); cOHdebug = cOH(:,:,i+1); cSO4debug = cSO4(:,:,i+1);
    cCadebug = cCa(:,:,i+1); cCaHCO3debug = cCaHCO3(:,:,i+1); cHCO3debug = cHCO3(:,:,i+1);
    cCO3debug = cCO3(:,:,i+1); 
    
    Idebug = 0.5*(2*cNaCl0 + cHSO4(:,:,i+1) + cH(:,:,i+1) + cOH(:,:,i+1) + 4*cSO4(:,:,i+1) + 4*cCa(:,:,i+1) +...
        cCaHCO3(:,:,i+1) + cHCO3(:,:,i+1) + 4*cCO3(:,:,i+1));
    I(:,:,i+1) = Idebug;
    
end


%% Sink/Source
SH = (cH(:,:,i) - cH0) / dt0;
SH2SO4 = (cH2SO4(:,:,i) - cH2SO40) / dt0;
SHSO4 = (cHSO4(:,:,i) - cHSO40) / dt0;
SSO4 = (cSO4(:,:,i) - cSO40) / dt0;
SCa = (cCa(:,:,i) - cCa0) / dt0;
SCaHCO3  = (cCaHCO3(:,:,i) - cCaHCO30) / dt0;
SCaT = (cCaT(:,:,i) - cCaT0) / dt0;
SH2CO3 = (cH2CO3(:,:,i) - cH2CO30) / dt0;
SHCO3 = (cHCO3(:,:,i) - cHCO30) / dt0;
SCO3 = (cCO3(:,:,i) - cCO30) / dt0;
SCT = (cCT(:,:,i) - cCT0) / dt0;

end
