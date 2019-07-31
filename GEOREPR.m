%% GEOthermal REinjection Lifetime PRediction (GEOREPR)

% This model (GEOREPR) is developed to predict geothermal reinjection lifetime limited 
% by the silica scaling under user-defined geothermal conditions at an acceptable 
% computational cost. 

% One can refer to the user guide to learn how to apply the model under different conditions. 

% This default example here models the decrement of injectivity due to silica scaling over 12 years   
%  in the Tiwi geothermal field, Philippines.                                                      

%% Start up
close all; clear all; clc; tic
GEOREPRStartUp();

%% Basic controlling inputs:
%     1. Initial conditions of the injectate (SGW)
%          and pre-existing fluid in the reservoir
[Tinj, pHinj, cDSinjppm, cPSinjppm] =...
    deal(160, 7, 1272, 0);
[Tres, pHres, cDSresppm, cPSresppm] =...
    deal(260, 7, 0, 0);
%     2.  Reinjection operating conditions
mDot = massflowrateSI(60, 'kg/s', Tinj);
%     3. Ageing
tAgeing = timeSI(35, 'min');
%     4. Timescale of interest
tmax = timeSI(12, 'a');
%     5. Reservoir geometry
[R0, H, a0mm, phiR0, Kf0] =...
    deal(0.122, 120, 7.7*2, 0.01, 5*10^-12);
%     6. Composition of injetate, initial reservoir fluids,
%          and the reservoir formation
[aHinj, cNaClinj] = deal(10^-pHinj, 0.09);
[cH2SO4inj, cHSO4inj, cSO4inj] = deal(0, 0, 0);
[cCainj, cCaCO30inj, cCaHCO3inj] = deal(0, 0, 0);
[cH2CO3inj, cHCO3inj, cCO3inj] = deal(0, 0, 0);
[aHres, cNaClres] = deal(10^-pHres, 0.07);
[cH2SO4res, cHSO4res, cSO4res] = deal(0, 0, 0);
[cCares, cCaCO30res, cCaHCO3res] = deal(0, 0, 0);
[cH2CO3res, cHCO3res, cCO3res] = deal(0, 0, 0);
[cHFinj, cFinj, cSiF6inj] = deal(0, 0, 0); [cHFres, cFres, cSiF6res] = deal(0, 0, 0);
[cHACinj, cACinj] = deal(0, 0); [cHACres, cACres] = deal(0, 0);
[sCaCO3, sInert] = deal(1, 0); 
%     8. Mesh quality
[Nx, Nz, Nt] = deal(50, 50, 50);

%% Input preprocessing 1/2
Ta = 0.5 * (Tinj + Tres); Taab = Ta + 273.15;
rhowater = -0.0025*Ta^2-0.1249*Ta+1005.2; % Density of water, kgm-3
xmax = 1.5 * (((mDot * tmax / rhowater + pi * R0^2 * H * phiR0) / (pi * H * phiR0)) ^ (1/2) - R0); 
zmax = lengthSI(a0mm/2, 'mm'); Rmax = 900; a0 = lengthSI(a0mm, 'mm');
zrockmax = (H - H * phiR0) / (H * phiR0 / lengthSI(a0mm, 'mm')) * 0.5;
rhorock = 2550; crock = 840; % Density and Specific heat of rock, kgm-3 & Jkg-1K-1, assumed
cwater = 0.0357*Ta^2-9.1868*Ta+4767.1; % Specific heat of water, Jkg-1K-1
Nzrock = 50; % Mesh quanlity of rock metrix
[x, dx, z, dz, zrock, dzrock, t0, ~] =...
    Mesh1(Nx, Nz, Nzrock, Nt, xmax, zmax, zrockmax, tmax);
dt = tmax/Nt; % time step 
t = 0; % Initial time, i.e. t0

%% Define the geometry and generate mesh
r = [R0 x'+R0 xmax+R0]; % cellfaces
zPS = [0 z]; z = [0 zrock zrockmax+z]; % cellfaces
m = createMeshCylindrical2D(r,z);
mPS = createMeshCylindrical2D(r,zPS);
% One-side fracture surface area
Abottom = 2*pi*m.cellcenters.x .* m.cellsize.x(2:end-1);
Abottom = [0; Abottom; 0];
% Control volumes
Ar = sCaCO3 .* [0; 2*pi.*m.cellcenters.x.*m.cellsize.x(2:end-1); 0];
Hr = m.cellcenters.y - m.cellcenters.y(Nzrock+1); Hr(1:Nzrock) = 0;
Hr = [0; Hr; 0]; Vr = Ar * Hr';

%% Constants
% DO NOT change the following constants at all times
g = 9.81; % m/s-2
kB = 1.38 * 10 ^ -23; %m2kgs-2K-1, Boltzmann constant
rhoSilica = 2196; %kg/m3, Density of amorphous silica
rhoNDS = 2.21 * 10^28; %m-3, The number density of SiO2 units in solids AS
vSIO2 = 0.00925*10^-3; %m3/mol, molar volume of silicon dioxide
Rg = 1.9872036; %cal/k/mol, gas constant
RSI = 8.3144598; %JK?1mol?1, gas constant
kBergs = 1.38 * 10^-16; %ergsK-1, Boltzmann constant
Fara = 96485.33289; %Cmol-1, Faraday constant
molarmassSilica = 0.06008; %kg mol-1
e = -1.602 * 10^-19; % C, Electron Charge
NA = 6.022 * 10^23; % Avogadro’s number
E = 8.8542E-12; % Fm-1, permittivity of free space

% MAY change the following constants if one is confident
Ew = relativepermittivityWater(Ta, 'degC'); % wiki
kwater = -5*10^-6*Ta^2+0.0013*Ta+0.6118; %Wm-1K-1, thermal conductivity of water
krock = 4.5; %Wm-1K-1, thermal conductivity of rock
mu = 5*10^-9*Ta^2-3*10^-6*Ta+0.0005; %Pas, Dynamic viscosity of injected and native reservoir fluids
nu = mu / rhowater; %m2s-1, Kinematic viscosity of injected and native reservoir fluids
lambdaNa0 = 50.11; %Scm2mol-1, conductivity of Na+ at 25 degC
lambdaH0 = 349.82; %Scm2mol-1, conductivity of H+
lambdaOH0 = 197.8; %Scm2mol-1, conductivity of OH-
lambdaCl0 = 76.34; %Scm2mol-1, conductivity of Cl- 
%http://www.currentseparations.com/issues/18-3/cs18-3c.pdf
lambdaSO40 = 160.0; %Scm2mol-1, conductivity of SO42-
lambdaAC0 = 40.9; %Scm2mol-1, conductivity of Ac-
% For most dilute natural
% waters, the conductivity of the complex can be
% ignored because the relative amount of complex formed is
% small and the mobility of the complex is low compared to
% that of the free ions. The protonated form of the sulfate
% anion (HSO4-) may also be safely ignored except in
% highly acidic waters containing sulfate, for example, acid
% mine drainage
lambdaHSO40 = 52; %Scm2mol-1, conductivity of HSO4-, http://kemia.ttk.pte.hu/pages/fizkem/oktatas/gyogyszeresz/CRC_Handbook_kivonat.pdf
lambdaCa0 = 119; %Scm2mol-1, conductivity of Ca2+
lambdaHCO30 = 44.5; %Scm2mol-1, conductivity of HCO3-, https://pubs.usgs.gov/wsp/2311/report.pdf
lambdaCO30 = 138.6; %Scm2mol-1, conductivity of Cl- 
lambdaCaHCO30 = 19; %Scm2mol-1, conductivity of Cl-, http://www.aqion.de/site/194
lambdaF0 = 55.4; %Scm2mol-1, conductivity of Cl- 
lambdaHF20 = 75; %Scm2mol-1, conductivity of Cl-, https://is.muni.cz/el/1431/podzim2016/C4020/um/pom/Ionic_Conductivity_and_Diffusion_at_Infinite_Dilution.pdf
lambdaSiF60 = 160; %Scm2mol-1, conductivity of Cl-, assumed
rSA = 4.565 * 10 ^ -10 / 2; %Scm2mol-1, conductivity of Cl- %averaged radius of Si(OH)4
% http://eprints.whiterose.ac.uk/119570/6/Supporting%20Information.pdf
rCaCO30 = (1.89 + 1.14) * 10^-10;
rH2CO3 = (2.02 + 0.3) * 10^-10; % assumed
rHF = 3.03 * 10^-10;
rHAC = 5.79 * 10^-10;
A = e^2 * (2*e^2*NA*rhowater/(E*Ew*kB*Taab))^0.5 / (2.302585*8*pi*E*Ew*kB*Taab); % factor of the Debye-Hückel equation
B = (2*e^2*NA*rhowater/(E*Ew*kB*Taab))^0.5 * 10^-9; % factor of the Debye-Hückel equation

%% Input preprocessing 2/2: check if equilibrium
% pseudo IS
ISinj =  0.5*(2*cNaClinj + cHSO4inj + 4*cSO4inj + 4*cCainj +...
    cCaHCO3inj + cHCO3inj + 4*cCO3inj + cFinj + 4*cSiF6inj);
ISres =  0.5*(2*cNaClres + cHSO4res + 4*cSO4res + 4*cCares +...
    cCaHCO3res + cHCO3res + 4*cCO3res + cFres + 4*cSiF6res);
Kw = 10 ^ -(6E-5*25^2-0.0286*25+14.463);
aOHinj = Kw / aHinj;
cOHinj = aOHinj / 10 ^ (-A * 1^2 * ISinj^0.5 / (1 + B * 3.5 * ISinj^0.5));
cHinj = aHinj / 10 ^ (-A * 1^2 * ISinj^0.5 / (1 + B * 9.0 * ISinj^0.5));
aOHres = Kw / aHres;
cOHres = aOHres / 10 ^ (-A * 1^2 * ISres^0.5 / (1 + B * 3.5 * ISres^0.5));
cHres = aHres / 10 ^ (-A * 1^2 * ISres^0.5 / (1 + B * 9.0 * ISres^0.5));
cSTinj = cH2SO4inj + cHSO4inj + cSO4inj;
cSTres = cH2SO4res + cHSO4res + cSO4res;
% H2SO4 related
% if cH2SO4inj == 0
%     [cHSO4inj, cSO4inj] =...
%         eq1H2SO4(Tinj, 'degC', cHinj, cSTinj, ISinj);
% else
%     [cH2SO4inj, cHinj, cHSO4inj, cSO4inj] =...
%         eqH2SO4(Tinj, 'degC', cHinj, cH2SO4inj, cHSO4inj, cSO4inj, ISinj);
% end
% [cHSO4res, cSO4res] =...
%     eq1H2SO4(Tres, 'degC', cHres, cSTres, ISres);
% HAC related
% [cHACinj, cHinj, cACinj] = eqHAC1(Tinj, 'degC', cHACinj, cHinj, cACinj, ISinj);
% [cHACres, cHres, cACres] = eqHAC1(Tres, 'degC', cHACres, cHres, cACres, ISres);
% HF related
% [cHFinj, cHinj, cFinj] = eqHF1(Tinj, 'degC', cHFinj, cHinj, cFinj, ISinj);
% [cHFres, cHres, cFres] = eqHF1(Tinj, 'degC', cHFres, cHres, cFres, ISres); %%
% update IS
ISinj =  0.5*(2*cNaClinj + cHSO4inj + cHinj + cOHinj + 4*cSO4inj + 4*cCainj +...
    cCaHCO3inj + cHCO3inj + 4*cCO3inj + cFinj + 4*cSiF6inj);
ISres =  0.5*(2*cNaClres + cHSO4res + cHres + cOHres + 4*cSO4res + 4*cCares +...
    cCaHCO3res + cHCO3res + 4*cCO3res + cFres + 4*cSiF6res);
cDSinj = cDSinjppm/1000/1000*rhowater; %kg/m3, ppm to kg/m3
cDSres = cDSresppm/1000/1000*rhowater; %kg/m3
cPSinj = cPSinjppm/1000/1000*rhowater; %kg/m3
cPSres = cPSresppm/1000/1000*rhowater; %kg/m3
% Ageing
if tAgeing == 0
    cDSw = cDSinj;
    cPSw = cPSinj;
    Astotali = zeros(length(r)+1,length(z)+1);
    Astotal0 = 0;
else
    [cDSw, Astotal0] = computePolymerisation(cDSinj, Tinj, 'degC', pHinj, ISinj, tAgeing, 0);
    % 0 at the end of computePolymerisation() means the esimation of
    % induction time is overridden
    cPSw = cDSinj - cDSw;
    Astotali = Astotal0 * ones(length(r)+1,length(z)+1);
end

%% Define boundary conditions
BCT = createBC(m); BCNaCl = createBC(m); BCH = createBC(m);
% BCOH = createBC(m);
BCDS = createBC(m); BCPS = createBC(m); BCAC = createBC(m);
BCHAC = createBC(m); BCCa = createBC(m); BCCaHCO3 = createBC(m);
BCCaCO30 = createBC(m); BCH2CO3 = createBC(m); BCHCO3 = createBC(m);
BCCO3 = createBC(m); BCHF = createBC(m); BCF = createBC(m); BCSiF6 = createBC(m);

% Define left boundary, others are 2nd type BC
BCT.left.a(:) = 0; BCT.left.b(:)=1; BCT.left.c(:) = Tinj;
BCNaCl.left.a(:) = 0; BCNaCl.left.b(:)=1; BCNaCl.left.c(:) = cNaClinj;
BCH.left.a(:) = 0; BCH.left.b(:)=1; BCH.left.c(:) = cHinj;
% BCOH.left.a(:) = 0; BCOH.left.b(:)=1; BCOH.left.c(:) = cOHinj;
BCDS.left.a(:) = 0; BCDS.left.b(:)=1; BCDS.left.c(:) = cDSw;
BCPS.left.a(:) = 0; BCPS.left.b(:)=1; BCPS.left.c(:) = cPSw;
BCHAC.left.a(:) = 0; BCHAC.left.b(:)=1; BCHAC.left.c(:) = cHACinj;
BCAC.left.a(:) = 0; BCAC.left.b(:)=1; BCAC.left.c(:) = cACinj;
BCCa.left.a(:) = 0; BCCa.left.b(:)=1; BCCa.left.c(:) = cCainj;
BCCaHCO3.left.a(:) = 0; BCCaHCO3.left.b(:)=1; BCCaHCO3.left.c(:) = cCaHCO3inj;
BCCaCO30.left.a(:) = 0; BCCaCO30.left.b(:)=1; BCCaCO30.left.c(:) = cCaCO30inj;
BCH2CO3.left.a(:) = 0; BCH2CO3.left.b(:)=1; BCH2CO3.left.c(:) = cH2CO3inj;
BCHCO3.left.a(:) = 0; BCHCO3.left.b(:)=1; BCHCO3.left.c(:) = cHCO3inj;
BCCO3.left.a(:) = 0; BCCO3.left.b(:)=1; BCCO3.left.c(:) = cCO3inj;
BCHF.left.a(:) = 0; BCHF.left.b(:)=1; BCHF.left.c(:) = cHFinj;
BCF.left.a(:) = 0; BCF.left.b(:)=1; BCF.left.c(:) = cFinj;
BCSiF6.left.a(:) = 0; BCSiF6.left.b(:)=1; BCSiF6.left.c(:) = cSiF6inj;

%% Define the transfer coeffs
% Diffusion field
% for heat transfer
DT = createCellVariable(m, 0);
DTf = kwater * ones(length(r) + 1, Nz + 1 - 1);
DTrock = krock * ones(length(r) + 1, Nzrock + 1 + 1);
DT.value = [DTrock DTf];
DTave = harmonicMean(DT);
MTdiff = diffusionTerm(DTave);
% for mass transfer
% in rock, D == 0
Dmassrock = zeros(length(r) + 1, Nzrock + 2);
% ions in fluids
Dionf = RSI / Fara^2 * ones(length(r) + 1, Nz);
% dissolved silica in fluids
DDSf = kB / (3 * pi * 2 * rSA * mu) * ones(length(r) + 1, Nz); % 10^-9 * has been deleted
% CaCO30
DCaCO30f = kB / (3 * pi * 2 * rCaCO30 * mu) * ones(length(r) + 1, Nz);
% H2CO3
DH2CO3f = kB / (3 * pi * 2 * rH2CO3 * mu) * ones(length(r) + 1, Nz);
% HF https://www.sciencedirect.com/science/article/pii/0013468685801791
DHFf = kB / (3 * pi * 2 * rHF * mu) * ones(length(r) + 1, Nz);
% HAC https://docplayer.net/14958547-Diffusion-and-fluid-flow.html
DHACf = kB / (3 * pi * 2 * rHAC * mu) * ones(length(r) + 1, Nz);

% Transient term coefficient
alfaT = createCellVariable(m, 1);
alfaTf = rhowater * cwater * ones(length(r) + 1, Nz + 1);
alfaTrock = rhorock * crock * ones(length(r) + 1, Nzrock + 1);
alfaT.value = [alfaTrock alfaTf];
alfamass = createCellVariable(m, 1);

% Velocity field
qDot = mDot / rhowater;
uTface = createFaceVariable(m, 0);
umassface = createFaceVariable(m, 0);
if m.dimension == 2
    tinj0 = m.cellcenters.x / (qDot / (H * phiR0 * w));
    umTf = rhowater * cwater * qDot / (H * phiR0 * w) * ones(size(m.facecenters.x));
elseif m.dimension == 2.5
    tinj0 = pi * H * phiR0 * (m.cellcenters.x.^2 - R0^2*ones(size(m.cellcenters.x))) / qDot;
    % Note: F = rho*cwater*u for heat transfer
    umTf = rhowater * cwater * mDot ./ m.facecenters.x / (2 * pi * H * rhowater * phiR0 * 1);
    % Note: F = u for mass transfer
end
ummassf = umTf / (rhowater * cwater);
h = lengthSI(a0mm/2, 'mm');
yf = 3/2*(1-(-m.cellcenters.y(Nzrock+2:end)+zrockmax+h).^2/h^2);
% yf = [zeros(length(z)-1-length(yf),1); ones(size(yf))]';
yf = [zeros(length(z)-1-length(yf),1); yf]';
uTface.xvalue = umTf * yf;
umassface.xvalue = ummassf * yf;
uPSface = createFaceVariable(m, 0);
uPSface.xvalue = umassface.xvalue;
MTconv =  convectionUpwindTerm(uTface); % Upwind is used for stability
% Mmassconv =  convectionTerm(umassface);
Mmassconv =  convectionUpwindTerm(umassface);

%% Initial values
% Heat
T0 = Tres;
Told = createCellVariable(m, T0, BCT);
T = Told;
[MTbc, RHSTbc] = boundaryCondition(BCT);
% Na+
cNaCl0 = cNaClres;
cNaClold = createCellVariable(m, cNaCl0, BCNaCl);
cNaCl = cNaClold;
[MNaClbc, RHSNaClbc] = boundaryCondition(BCNaCl);
% % HSO4
% cH2SO4old = createCellVariable(m, 0);
% cHSO40 = cHSO4res;
% cHSO4old = createCellVariable(m, cHSO40, BCHAC);
% cHSO4 = cHSO4old;
% [MHSO4bc, RHSHSO4bc] = boundaryCondition(BCHAC);
% HAC
% cH2SO4old = createCellVariable(m, 0);
cHAC0 = cHACres;
cHACold = createCellVariable(m, cHAC0, BCHAC);
cHAC = cHACold;
[MHACbc, RHSHACbc] = boundaryCondition(BCHAC);
% % SO4
% cSO40 = cSO4res;
% cSO4old = createCellVariable(m, cSO40, BCAC);
% cSO4 = cSO4old;
% [MSO4bc, RHSSO4bc] = boundaryCondition(BCAC);
cAC0 = cACres;
cACold = createCellVariable(m, cAC0, BCAC);
cAC = cACold;
[MACbc, RHSACbc] = boundaryCondition(BCAC);
% Ca
cCa0 = cCares;
cCaold = createCellVariable(m, cCa0, BCCa);
cCa = cCaold;
[MCabc, RHSCabc] = boundaryCondition(BCCa);
% CaHCO3
cCaHCO30 = cCaHCO3res;
cCaHCO3old = createCellVariable(m, cCaHCO30, BCCaHCO3);
cCaHCO3 = cCaHCO3old;
[MCaHCO3bc, RHSCaHCO3bc] = boundaryCondition(BCCaHCO3);
% CaT
cCaTold = createCellVariable(m,0);
% H2CO3
cH2CO30 = cH2CO3res;
cH2CO3old = createCellVariable(m, cH2CO30, BCH2CO3);
cH2CO3 = cH2CO3old;
[MH2CO3bc, RHSH2CO3bc] = boundaryCondition(BCH2CO3);
% HCO3
cHCO30 = cHCO3res;
cHCO3old = createCellVariable(m, cHCO30, BCHCO3);
cHCO3 = cHCO3old;
[MHCO3bc, RHSHCO3bc] = boundaryCondition(BCHCO3);
% CO3
cCO30 = cCO3res;
cCO3old = createCellVariable(m, cCO30, BCCO3);
cCO3 = cCO3old;
[MCO3bc, RHSCO3bc] = boundaryCondition(BCCO3);
% CT
cCTold = createCellVariable(m,0);
% H+
cH0 = cHres;
cHold = createCellVariable(m, cH0, BCH);
aHold = createCellVariable(m, 0);
pHold = createCellVariable(m, 0);
cH = cHold;
[MHbc, RHSHbc] = boundaryCondition(BCH);
% OH-
cOHold = createCellVariable(m, 0);
aOHold = createCellVariable(m, 0);
% cOH0 = cOHres;
% cOHold = createCellVariable(m, cOH0, BCOH);
% cOH = cOHold;
% [MOHbc, RHSOHbc] = boundaryCondition(BCOH);
% DS
cDS0 = cDSres;
cDSold = createCellVariable(m, cDS0, BCDS);
cDSold.value(:,1:Nzrock+1) = 0;
cDS = cDSold;
[MDSbc, RHSDSbc] = boundaryCondition(BCDS);
% PS
cPS0 = cPSres;
cPSold = createCellVariable(m, cPS0, BCPS);
cPS = cPSold;
[MPSbc, RHSPSbc] = boundaryCondition(BCPS);
% HF
cHF0 = cHFres;
cHFold = createCellVariable(m, cHF0, BCHF);
cHF = cHFold;
[MHFbc, RHSHFbc] = boundaryCondition(BCHF);
% F-
cF0 = cFres;
cFold = createCellVariable(m, cF0, BCF);
cF = cFold;
[MFbc, RHSFbc] = boundaryCondition(BCF);
% SiF62-
cSiF60 = cSiF6res;
cSiF6old = createCellVariable(m, cSiF60, BCSiF6);
cSiF6 = cSiF6old;
[MSiF6bc, RHSSiF6bc] = boundaryCondition(BCSiF6);
% IS
ISold = createCellVariable(m, ISres);
ISpos = ISold.value; ISpos(ISpos<0) = 0;
% jDep
jDep = createCellVariable(m, 0);
% TVD scheme
FL = fluxLimiter('Superbee');

%% solver
figure(1);
figtitle = sgtitle('GEOREPR Simulation results t = 0 s');
hdep = zeros(size(m.cellcenters.x));
i = 1;

while t < tmax
    disp(t)
    t = t + dt;
    tinj = tinj0 + tAgeing;
    if t < 60
        dispt = t;
        dispunit = 'seconds';
    elseif t < 3600
        dispt = t/60;
        dispunit = 'minutes';
    elseif t < 3600*24
        dispt = t/3600;
        dispunit = 'hours';
    elseif t < 3600*24*30
        dispt = t/3600/24;
        dispunit = 'days';
    elseif t < 3600*24*30*12
        dispt = t/3600/24/30;
        dispunit = 'months';
    else
        dispt = t/3600/24/365;
        dispunit = 'years';
    end
    %% Heat transfer
    [MTtrans, RHSTtrans] = transientTerm(Told, dt, alfaT);
    MT = MTtrans-MTdiff+MTconv+MTbc;
    RHST = RHSTtrans+RHSTbc;
    T = solvePDE(m, MT, RHST);
    Told = T;
    Tabold = Told.value + 273.15;
    %% Background eletrolytes
    %% Na+
    lambdaNa = (1.258 * Told.value - 35.808 + (0.0044 * Told.value +...
        0.8841) * lambdaNa0) * 10^-4;
    DNaCl = createCellVariable(m, 0);
    DNaCl.value = [Dmassrock Dionf] .* Tabold .* lambdaNa;
    DNaClave = harmonicMean(DNaCl);
    MNaCldiff = diffusionTerm(DNaClave);
    [MNaCltrans, RHSNaCltrans] = transientTerm(cNaClold, dt, alfamass);
    %% H2SO4, SC for sulphuric acid, cSC is always zero
    %% HAC
%     lambdaHSO4 = (1.258 * Told.value - 35.808 + (0.0044 * Told.value +...
%         0.8841) * lambdaHSO40) * 10^-4;
    DHAC = createCellVariable(m, 0);
%     DHAC.value = [Dmassrock Dionf] .* Tabold .* lambdaHSO4;
    DHAC.value = [Dmassrock DHACf] .* Tabold;
    DHACave = harmonicMean(DHAC);
    MHACdiff = diffusionTerm(DHACave);
    [MHSO4trans, RHSHSO4trans] = transientTerm(cHACold, dt, alfamass);
    %% AC
    lambdaAC = (1.258 * Told.value - 35.808 + (0.0044 * Told.value +...
        0.8841) * lambdaAC0) * 10^-4;
    DAC = createCellVariable(m, 0);
    DAC.value = [Dmassrock Dionf] .* Tabold .* lambdaAC;
    DACave = harmonicMean(DAC);
    MACdiff = diffusionTerm(DACave);
    [MACtrans, RHSACtrans] = transientTerm(cACold, dt, alfamass);
    %% Ca2+
    lambdaCa = (1.258 * Told.value - 35.808 + (0.0044 * Told.value +...
        0.8841) * lambdaCa0) * 10^-4;
    DCa = createCellVariable(m, 0);
    DCa.value = [Dmassrock Dionf] .* Tabold .* lambdaCa;
    DCaave = harmonicMean(DCa);
    MCadiff = diffusionTerm(DCaave);
    [MCatrans, RHSCatrans] = transientTerm(cCaold, dt, alfamass);
    MCa = MCatrans-MCadiff+Mmassconv+MCabc;
    RHSCa = RHSCatrans+RHSCabc;
    cCa = solvePDE(m, MCa, RHSCa);
    cCaold = cCa;
    %% CaHCO3+
    lambdaCaHCO3 = (1.258 * Told.value - 35.808 + (0.0044 * Told.value +...
        0.8841) * lambdaCaHCO30) * 10^-4;
    DCaHCO3 = createCellVariable(m, 0);
    DCaHCO3.value = [Dmassrock Dionf] .* Tabold .* lambdaCaHCO3;
    DCaHCO3ave = harmonicMean(DCaHCO3);
    MCaHCO3diff = diffusionTerm(DCaHCO3ave);
    [MCaHCO3trans, RHSCaHCO3trans] = transientTerm(cCaHCO3old, dt, alfamass);
    %% H2CO3
    DH2CO3 = createCellVariable(m, 0);
    DH2CO3.value = [Dmassrock DH2CO3f] .* Tabold;
    DH2CO3ave = harmonicMean(DH2CO3);
    MH2CO3diff = diffusionTerm(DH2CO3ave);
    [MH2CO3trans, RHSH2CO3trans] = transientTerm(cH2CO3old, dt, alfamass);
    %% HCO3-
    lambdaHCO3 = (1.258 * Told.value - 35.808 + (0.0044 * Told.value +...
        0.8841) * lambdaHCO30) * 10^-4;
    DHCO3 = createCellVariable(m, 0);
    DHCO3.value = [Dmassrock Dionf] .* Tabold .* lambdaHCO3;
    DHCO3ave = harmonicMean(DHCO3);
    MHCO3diff = diffusionTerm(DHCO3ave);
    [MHCO3trans, RHSHCO3trans] = transientTerm(cHCO3old, dt, alfamass);
    %% CO32-
    lambdaCO3 = (1.258 * Told.value - 35.808 + (0.0044 * Told.value +...
        0.8841) * lambdaCO30) * 10^-4;
    DCO3 = createCellVariable(m, 0);
    DCO3.value = [Dmassrock Dionf] .* Tabold .* lambdaCO3;
    DCO3ave = harmonicMean(DCO3);
    MCO3diff = diffusionTerm(DCO3ave);
    [MCO3trans, RHSCO3trans] = transientTerm(cCO3old, dt, alfamass);
    %% H+
    lambdaH = (1.258 * Told.value - 35.808 + (0.0044 * Told.value + 0.8841) *...
        lambdaH0) * 10^-4;
    DH = createCellVariable(m, 0);
    DH.value = [Dmassrock Dionf] .* Tabold .* lambdaH;
    DHave = harmonicMean(DH);
    MHdiff = diffusionTerm(DHave);
    [MHtrans, RHSHtrans] = transientTerm(cHold, dt, alfamass);
    %% OH-
    %     lambdaOH = (1.258 * Told.value - 35.808 + (0.0044 * Told.value + 0.8841) *...
    %         lambdaOH0) * 10^-4;
    %     DOH = createCellVariable(m, 0);
    %     DOH.value = [Dmassrock Dionf] .* Tabold .* lambdaOH;
    %     DOHave = harmonicMean(DOH);
    %     MOHdiff = diffusionTerm(DOHave);
    %     [MOHtrans, RHSOHtrans] = transientTerm(cOHold, dt, alfamass);
    %% HF
    DHF = createCellVariable(m, 0);
    DHF.value = [Dmassrock DHFf] .* Tabold;
    DHFave = harmonicMean(DHF);
    MHFdiff = diffusionTerm(DHFave);
    [MHFtrans, RHSHFtrans] = transientTerm(cHFold, dt, alfamass);
    %% F-
    lambdaF = (1.258 * Told.value - 35.808 + (0.0044 * Told.value + 0.8841) *...
        lambdaF0) * 10^-4;
    DF = createCellVariable(m, 0);
    DF.value = [Dmassrock Dionf] .* Tabold .* lambdaF;
    DFave = harmonicMean(DF);
    MFdiff = diffusionTerm(DFave);
    [MFtrans, RHSFtrans] = transientTerm(cFold, dt, alfamass);
    %% SiF62-
    lambdaSiF6 = (1.258 * Told.value - 35.808 + (0.0044 * Told.value + 0.8841) *...
        lambdaSiF60) * 10^-4;
    DSiF6 = createCellVariable(m, 0);
    DSiF6.value = [Dmassrock Dionf] .* Tabold .* lambdaSiF6;
    DSiF6ave = harmonicMean(DSiF6);
    MSiF6diff = diffusionTerm(DSiF6ave);
    [MSiF6trans, RHSSiF6trans] = transientTerm(cSiF6old, dt, alfamass);
    %% Combined geochem subsolver
    SH = createCellVariable(m, 0); 
%     SHAC = createCellVariable(m, 0);
%     SAC = createCellVariable(m, 0); SCa = createCellVariable(m, 0);
%     SCaHCO3 = createCellVariable(m, 0); SH2CO3 = createCellVariable(m, 0); 
%     SHCO3 = createCellVariable(m, 0); SCO3 = createCellVariable(m, 0);
%     SHF0 = createCellVariable(m, 0); SF = createCellVariable(m, 0); 
%     SSiF6 = createCellVariable(m, 0); 
%     
% %     [qH, qH2SO4, qHSO4, qSO4, qCa, qCaHCO3, qH2CO3, qHCO3, qCO3] =...
% %         SCalcite1(Told.value, 'degC',...
% %         cNaClold.value, cHold.value, cH2SO4old.value, cHSO4old.value,...
% %         cSO4old.value, cCaold.value,cCaHCO3old.value, cCaTold.value,...
% %         cH2CO3old.value, cHCO3old.value, cCO3old.value, cCTold.value, m, a0, sCaCO3, dt);
%     
%     [qH, qHAC, qAC, qCa, qCaHCO3, qH2CO3, qHCO3, qCO3] = SCalciteAC(Told.value, 'degC',...
%         cNaClold.value, cHold.value, cHACold.value, cACold.value, cCaold.value,...
%         cCaHCO3old.value, cCaTold.value, cH2CO3old.value, cHCO3old.value,...
%         cCO3old.value, cCTold.value, m, a0, sCaCO3, dt);
% 
%     qH(1,:) = 0; qH(:,1) = 0; qH(end,:) = 0; qH(:,end) = 0; qH(isnan(qH)) = 0;
% %     qHSO4(1,:) = 0; qHSO4(:,1) = 0; qHSO4(end,:) = 0; qHSO4(:,end) = 0; qHSO4(isnan(qHSO4)) = 0;
% %     qSO4(1,:) = 0; qSO4(:,1) = 0; qSO4(end,:) = 0; qSO4(:,end) = 0; qSO4(isnan(qSO4)) = 0;
%     qHAC(1,:) = 0; qHAC(:,1) = 0; qHAC(end,:) = 0; qHAC(:,end) = 0; qHAC(isnan(qHAC)) = 0;
%     qAC(1,:) = 0; qAC(:,1) = 0; qAC(end,:) = 0; qAC(:,end) = 0; qAC(isnan(qAC)) = 0;
%     qCa(1,:) = 0; qCa(:,1) = 0; qCa(end,:) = 0; qCa(:,end) = 0;
%     qCaHCO3(1,:) = 0; qCaHCO3(:,1) = 0; qCaHCO3(end,:) = 0;
%     qCaHCO3(:,end) = 0; qCaHCO3(isnan(qCaHCO3)) = 0;
%     qH2CO3(1,:) = 0; qH2CO3(:,1) = 0; qH2CO3(end,:) = 0; qH2CO3(:,end) = 0; qH2CO3(isnan(qH2CO3)) = 0;
%     qHCO3(1,:) = 0; qHCO3(:,1) = 0; qHCO3(end,:) = 0; qHCO3(:,end) = 0; qHCO3(isnan(qHCO3)) = 0;
%     qCO3(1,:) = 0; qCO3(:,1) = 0; qCO3(end,:) = 0; qCO3(:,end) = 0; qCO3(isnan(qCO3)) = 0;
%     
%     qH = -qH; qHAC = -qHAC; qAC = -qAC; qCa = -qCa; qCaHCO3 = -qCaHCO3; 
%     qH2CO3 = -qH2CO3; qHCO3 = -qHCO3; qCO3 = -qCO3; 
% %     qH(:,1:Nzrock+2) = 0; qH(:,Nzrock+4:end) = 0;
%     
%     [qHF, qH1, qF, qSiF6] = SHF1(Told.value, 'degC', cHFold.value, cHold.value, cFold.value,...
%         cSiF6old.value, ISold.value, m, a0, dt);
%     
%     qHF(1,:) = 0; qHF(:,1) = 0; qHF(end,:) = 0; qHF(:,end) = 0; qHF(isnan(qHF)) = 0;
%     qH1(1,:) = 0; qH1(:,1) = 0; qH1(end,:) = 0; qH1(:,end) = 0; qH1(isnan(qH1)) = 0;
%     qF(1,:) = 0; qF(:,1) = 0; qF(end,:) = 0; qF(:,end) = 0; qF(isnan(qF)) = 0;
%     qSiF6(1,:) = 0; qSiF6(:,1) = 0; qSiF6(end,:) = 0; qSiF6(:,end) = 0; qSiF6(isnan(qSiF6)) = 0;
% %     qH = qH - qH1; 
%     qHF = -qHF; qF = -qF; qSiF6 = -qSiF6;
% 
%     SH.value = qH; SHAC.value = qHAC; SAC.value = qAC; 
%     SCa.value = qCa; SCaHCO3.value = qCaHCO3;
%     SH2CO3.value = qH2CO3; SHCO3.value = qHCO3; SCO3.value = qCO3;
%     
% %     SH.value = -qH1;
%     SHF0.value = qHF; SF.value = qF; SSiF6.value = qSiF6;
    MHlin = linearSourceTerm(SH); RHSHlin = constantSourceTerm(SH);
%     MHAClin = linearSourceTerm(SHAC); RHSHAClin = constantSourceTerm(SHAC);
%     MAClin = linearSourceTerm(SAC); RHSAClin = constantSourceTerm(SAC);
%     MCalin = linearSourceTerm(SCa); RHSCalin = constantSourceTerm(SCa);
%     MCaHCO3lin = linearSourceTerm(SCaHCO3); RHSCaHCO3lin = constantSourceTerm(SCaHCO3);
%     MH2CO3lin = linearSourceTerm(SH2CO3); RHSH2CO3lin = constantSourceTerm(SH2CO3);
%     MHCO3lin = linearSourceTerm(SHCO3); RHSHCO3lin = constantSourceTerm(SHCO3);
%     MCO3lin = linearSourceTerm(SCO3); RHSCO3lin = constantSourceTerm(SCO3);
%     MHFlin = linearSourceTerm(SHF0); RHSHFlin = constantSourceTerm(SHF0);
%     MFlin = linearSourceTerm(SF); RHSFlin = constantSourceTerm(SF);
%     MSiF6lin = linearSourceTerm(SSiF6); RHSSiF6lin = constantSourceTerm(SSiF6);
%   
%     % NaCl
%     MNaCl = MNaCltrans-MNaCldiff+Mmassconv+MNaClbc;
%     RHSNaCl = RHSNaCltrans+RHSNaClbc;
%     cNaCl = solvePDE(m, MNaCl, RHSNaCl);
%     cNaClold = cNaCl;
    % H
    MH = MHtrans-MHdiff+Mmassconv+MHbc + MHlin;
    RHSH = RHSHtrans+RHSHbc - RHSHlin;
    cH = solvePDE(m, MH, RHSH);
    cHold = cH;
%     
%     % HAC
%     MHAC = MHSO4trans-MHACdiff+Mmassconv+MHACbc + MHAClin;
%     RHSHAC = RHSHSO4trans+RHSHACbc - RHSHAClin;
%     cHAC = solvePDE(m, MHAC, RHSHAC);
%     cHACold = cHAC;
%     % AC
%     MAC = MACtrans-MACdiff+Mmassconv+MACbc + MAClin;
%     RHSAC = RHSACtrans+RHSACbc - RHSAClin;
%     cAC = solvePDE(m, MAC, RHSAC);
%     cACold = cAC;
%     
% %     % HAC
% %     MHAC = MHACtrans-MHACdiff+Mmassconv+MHACbc + MHAClin;
% %     RHSHAC = RHSHACtrans+RHSHACbc - RHSHAClin;
% %     cHAC = solvePDE(m, MHAC, RHSHAC);
% %     cHACold = cHAC;
% %     
% %     %AC
% %     MAC = MACtrans-MACdiff+Mmassconv+MACbc + MAClin;
% %     RHSAC = RHSACtrans+RHSACbc - RHSAClin;
% %     cAC = solvePDE(m, MAC, RHSAC);
% %     cACold = cAC;
%     
%     % Ca ???
%     MCa = MCatrans-MCadiff+Mmassconv+MCabc + MCalin;
%     RHSCa = RHSCatrans+RHSCabc - RHSCalin;
%     cCa = solvePDE(m, MCa, RHSCa);
%     cCaold = cCa;
%     % CaHCO3
%     MCaHCO3 = MCaHCO3trans-MCaHCO3diff+Mmassconv+MCaHCO3bc + MCaHCO3lin;
%     RHSCaHCO3 = RHSCaHCO3trans+RHSCaHCO3bc - RHSCaHCO3lin;
%     cCaHCO3 = solvePDE(m, MCaHCO3, RHSCaHCO3);
%     cCaHCO3old = cCaHCO3;
%     % H2CO3
%     MH2CO3 = MH2CO3trans-MH2CO3diff+Mmassconv+MH2CO3bc + MH2CO3lin;
%     RHSH2CO3 = RHSH2CO3trans+RHSH2CO3bc - RHSH2CO3lin;
%     cH2CO3 = solvePDE(m, MH2CO3, RHSH2CO3);
%     cH2CO3old = cH2CO3;
%     % HCO3
%     MHCO3 = MHCO3trans-MHCO3diff+Mmassconv+MHCO3bc + MHCO3lin;
%     RHSHCO3 = RHSHCO3trans+RHSHCO3bc - RHSHCO3lin;
%     cHCO3 = solvePDE(m, MHCO3, RHSHCO3);
%     cHCO3old = cHCO3;
%     % CO32
%     MCO3 = MCO3trans-MCO3diff+Mmassconv+MCO3bc + MCO3lin;
%     RHSCO3 = RHSCO3trans+RHSCO3bc - RHSCO3lin;
%     cCO3 = solvePDE(m, MCO3, RHSCO3);
%     cCO3old = cCO3;
%     % CT
%     cCTold = cH2CO3old + cHCO3old + cCO3old + cCaHCO3old;
%     % Ca total
%     cCaTold = cCaold + cCaHCO3old;
%     % HF
%     MHF = MHFtrans-MHFdiff+Mmassconv+MHFbc + MHFlin;
%     RHSHF = RHSHFtrans+RHSHFbc - RHSHFlin;
%     cHF = solvePDE(m, MHF, RHSHF);
%     cHFold = cHF;
%     % F
%     MF = MFtrans-MFdiff+Mmassconv+MFbc + MFlin;
%     RHSF = RHSFtrans+RHSFbc - RHSFlin;
%     cF = solvePDE(m, MF, RHSF);
%     cFold = cF;
%     % SiF6
%     MSiF6 = MSiF6trans-MSiF6diff+Mmassconv+MSiF6bc + MSiF6lin;
%     RHSSiF6 = RHSSiF6trans+RHSSiF6bc - RHSSiF6lin;
%     cSiF6 = solvePDE(m, MSiF6, RHSSiF6);
%     cSiF6old = cSiF6;
%     
%     [cHFold.value, cHold.value, cFold.value] = eqHF1(Told.value, 'degC', cHFold.value,...
%         cHold.value, cFold.value, ISold.value);
%     
% %     [cHACold.value, cHold.value, cACold.value] = eqHAC1(Told.value, 'degC', cHACold.value,...
% %         cHold.value, cACold.value, ISold.value);
%     
    gammaHold = 10 .^ (-A * 1^2 * ISpos.^0.5 ./ (1 + B * 9.0 * ISpos.^0.5));
%     
    aHold.value = cHold.value .* gammaHold;
%     aOHold.value = Kw ./ aHold.value;
%     cOHold.value = aOHold.value / 10 ^ (-A * 1^2 * ISinj^0.5 / (1 + B * 3.5 * ISinj^0.5));
    pHold.value = - log10(aHold.value);
%     
%     % IS
%     ISold.value = 0.5*(2*cNaClold.value + cHACold.value + cHold.value +...
%         cOHold.value + 4*cACold.value + 4*cCaold.value +...
%         cCaHCO3old.value + cHCO3old.value + 4*cCO3old.value + cFold.value + 4*cSiF6old.value);
    % Override Chem Model
    ISold.value = 0.5*2*cNaClold.value;
    ISpos = ISold.value; ISpos(ISpos<0) = 0;
    %% DS
    DDS = createCellVariable(m, 0);
    DDS.value = [Dmassrock DDSf] .* Tabold;
    DDSave = harmonicMean(DDS);
    MDSdiff = diffusionTerm(DDSave);
    [MDStrans, RHSDStrans] = transientTerm(cDSold, dt, alfamass);
    % DS source
    SDS = createCellVariable(m, 0);
    if mDot == 0
        [qDS, Astotali] = sinkDS(cDSinj*[zeros(length(r) + 1, Nzrock + 2) ones(length(r) + 1, Nz)], cDSold.value,...
            Astotali, Told.value, 'degC', pHold.value, ISold.value, dt);
    else
        Astotali = asTotal1(cDSinj, cDSw, Astotal0, Told.value, 'degC', pHold.value, ISold.value, tinj);
        Astotali = [0 Astotali 0]' * [zeros(1, Nzrock + 2) ones(1, Nz)];
        [qDS, ~] = sinkDS(cDSinj*[zeros(length(r) + 1, Nzrock + 2) ones(length(r) + 1, Nz)], cDSold.value,...
            Astotali, Told.value, 'degC', pHold.value, ISold.value, dt);
    end
    qDS = -qDS;
    % Assign 0 sink to ghost cells
    qDS(1,:) = 0; qDS(:,1) = 0; qDS(end,:) = 0; qDS(:,end) = 0;
    Rmd = computeRmd2(cDSold.value, Told.value, 'degC', pHold.value, ISold.value);
    pHbottom = pHold.value(:, Nzrock+3); ISbottom = ISold.value(:, Nzrock+3); Tbottom = Told.value(:, Nzrock+3);
    qMD = 2*Rmd / a0; % As = 2/a0
    qMD(1,:) = 0; qMD(:,1) = 0; qMD(end,:) = 0; qMD(:,end) = 0;
    SDS.value = qDS + qMD;
    SDS.value(isnan(SDS.value)) = 0;
    cDSmin = solubilityAS(Told.value, 'degC', pHold.value, ISold.value);
    if t > 0
        ifSDS = ((cDSold.value - SDS.value*dt) >= cDSmin);
        SDS.value = ifSDS.*SDS.value + (ones(size(ifSDS))-ifSDS).*qDS;
    end
    MDSlin = linearSourceTerm(SDS);
    RHSDSlin = constantSourceTerm(SDS);
    MDS = MDStrans - MDSdiff + Mmassconv + MDSbc + MDSlin;
    RHSDS = RHSDStrans + RHSDSbc - RHSDSlin;
    cDS = solvePDE(m, MDS, RHSDS);
    cDSold = cDS;
    %% PS
    %Redefine velocity field for PS
    rPS = particleRadiusM(cDSinj, Told.value, 'degC', pHold.value, ISold.value, tinj) * 10^-9;
    DPSf = kB ./ (3 * pi * 2 * rPS * mu);
    DPSf(DPSf==inf) = 0;
    DPSf = DPSf(:,Nzrock + 3:end);
    DPS = createCellVariable(m, 0);
    DPS.value = [Dmassrock DPSf] .* Tabold; DPSave = harmonicMean(DPS); 
    DPSave.xvalue(isnan(DPSave.xvalue)) = 0;
    DPSave.yvalue(isnan(DPSave.yvalue)) = 0;
    MPSdiff = diffusionTerm(DPSave);
    TaboldF = Tabold;
    TaboldF(1,:) = []; TaboldF(end,:) = []; TaboldF(:,1)=[];
    F = zeros(size(TaboldF));
    uPSface.yvalue = F .* DPSave.yvalue ./ (kB * TaboldF);
    uPSface.yvalue(isnan(uPSface.yvalue)) = 0; %?
    MPSconv =  convectionTerm(uPSface);
    [MPStrans, RHSPStrans] = transientTerm(cPSold, dt, alfamass);
    %PS source
    SPS = createCellVariable(m, 0);
    cPSbottom = cPSold.value(:, Nzrock+3);
    DPSbottom = DPS.value(:, Nzrock+3);
    Tabbottom = Tabold(:, Nzrock+3);
    surfaceChargebottom = computeSurfacecharge(pHbottom, ISbottom);
    dnmbottom = 2*10^9*rPS(:, Nzrock+3);
    Wbottom = computeStability1(dnmbottom, pHbottom, ISbottom);
    qCD = -4*kB*Tabbottom ./ (3*mu*Wbottom.*Abottom) .*...
        (1 - mPS.cellsize.y(1)./(2*DPSbottom).*(4*kB*Tabbottom./(3*mu*Wbottom.*Abottom))).^(-2).*...
        (-4*kB*Tabbottom ./ (3*mu*Wbottom.*Abottom) ./ (2*DPSbottom)) .* cPSbottom .*...
        [zeros(1, Nzrock + 2) ones(1, 1) zeros(1, Nz-1)];
    qCD(isnan(qCD)) = 0;
    qCD(1,:) = 0; qCD(:,1) = 0; qCD(end,:) = 0; qCD(:,end) = 0;
    SPS.value = -qDS + qCD;
    MPSlin = linearSourceTerm(SPS);
    RHSPSlin = constantSourceTerm(SPS);
    MPS = MPStrans - MPSdiff + MPSconv + MPSbc + MPSlin;
    RHSPS = RHSPStrans + RHSPSbc - RHSPSlin;
    cPS = solvePDE(m, MPS, RHSPS);
    cPSold = cPS;
    %% Deposition
    jCDM = 4*kB*Tabbottom ./ (3*mu*Wbottom.*Abottom) .* cPSbottom;
    jCDM(1) = []; jCDM(end) = []; jCDM(isnan(jCDM)) = 0;
    jMDM = Rmd(:,Nzrock+3); jMDM(1) = []; jMDM(end) = [];
    jDepM = jCDM + jMDM;
    jDep.value = [[0; jDepM; 0] zeros(size(cDS.value(:,1:end-1)))];
    kDepM = jDepM/rhoSilica;
    hdep = hdep + 2*kDepM*dt;
    Vdep(i) = ones(size(m.cellcenters.x'))*(2*pi*(m.cellcenters.x).*m.cellsize.x(2:end-1).* hdep);
    MdisCa = rhowater * cCaTold.value .* Vr; MdisCa(isnan(MdisCa)) = 0;
    MdisSi = rhowater * cSiF6old.value .* Vr; MdisSi(isnan(MdisSi)) = 0;
    VdisCa(i) = sum(sum(MdisCa)) * 100.0869 / 1000 / 2711; % m3
    VdisSi(i) = sum(sum(MdisSi)) * 60.08 / 1000 / 2196; % m3
    xinj(i) = (((mDot * t / rhowater + pi * R0^2 * H * phiR0) / (pi * H * phiR0)) ^ (1/2) - R0);
    tyear(i) = t/3600/24/365;
    i = i + 1;
    %% Plot
    set(gcf, 'Position',  [150, 150, 750, 600]);
    subplot(4,1,1);
    visualizeCells(Told);
    colormap(gca, 'hot')
    shading interp
    xlabel('Radial distance (m)')
    ylabel('Vertial distance (m)')
    title1 = ['Temperature (' char(176) 'C)'];
    title(title1);
    subplot(4,1,2);
    visualizeCells(cDSold);
    shading interp
    xlabel('Radial distance (m)')
    ylabel('Vertial distance (m)')
    title('Concentration of dissolved silica (kg m^{-3})')
    subplot(4,1,3);
    visualizeCells(cPSold);
    colormap(gca, 'winter')
    shading interp
    xlabel('Radial distance (m)')
    ylabel('Vertial distance (m)')
    title('Concentration of polymerised silica (kg m^{-3})')
    subplot(4,1,4);
    visualizeCells(jDep);
    colormap(gca, flipud(gray))
    shading interp
    xlabel('Radial distance (m)')
    ylabel('Vertial distance (m)')
    title('Silica deposition rate (kg m^{-2} s^{-1})')
    figtitle.String = ['GEOREPR Simulation results t = ', sprintf('%d %s',dispt,dispunit)];
    drawnow
end

%% Post Processing
% Observation
tob = [0.2 0.3 2.1 2.6 2.9 4 4.4 5.3 7.25 10.8 10.85 11.45];
IIob = [4.45 4.2 2.55 2.45 1.4 1.6 1.35 1.5 0.5 0.3 0.25 0.05];

phif = ones(size(Vdep)) - Vdep / (pi*(Rmax^2-R0^2)*a0(1));
phic = 0.98; n =  4;
phif0 = ones(size(phif));
Kf = Kf0 * ((phif - phic) ./ (phif0 - phic)).^n;
II = 10^5 * rhowater * 2 * pi * Kf * a0 / (mu * log(Rmax / R0)) * H * phiR0 / a0;
IIi = smooth(II,10);

figure(2)
hold on
plot(tyear, IIi);
plot(tob, IIob, 's');
hold off
grid on
legend('Prediction','Observation');
xlabel('Time (year)','fontsize',12,'fontweight','light');
ylabel('Injectivity index (kgs^{-1}bar^{-1})','fontsize',12,'fontweight','light');

toc

