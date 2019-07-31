function W = computeStability1(dnm, pH0, IS0)

% clear all; close all; clc;
% dnm = 100;
% pH0 = 8;
% IS0 = 0.01;

[Nx, Ny] = size(dnm);

for m = 1:Nx
    for n = 1:Ny
        
        ctd0 = 1;
      %% Enter the inputs (classic DLVO theory)
        TdegC = 25; % degC, temperature
        IS = IS0(m,n); % M, ionic strength
        pH = pH0(m,n); %
        d1nm = dnm(m,n); % Diametre of particle 1
        d2nm = 1000; % Diametre of particle 2
        
       %% Enter the inputs (Soft particle model)
        
        % Hsminnm = 30.66; % Separation distance where aggregation happens
%         d0nm = 10.^(-(0.004075*pH + 0.5511) .* log10(IS) + (0.07238*pH - 0.6817));
        d0nm = computed0nm(d1nm, pH, IS);
%         d0nm = 7.63;
        Hsminnm = 2*d0nm(ctd0); % Capture distance ##
        sigma0mC = - computeSurfacecharge(pH, IS); % mC/m2, a function of pH and IS
        d0(ctd0) = d0nm(ctd0) * 10^-9; % Thickness of gel-layer
        sigma0 = sigma0mC * 10^-3;
        
      %% Constants
        kB = 1.38 * 10^-23; %m2kgs-2K-1, Boltzmann constant
        n1 = 1.43; % Refractive Index of Silica
        n3 = 1.33; % Refractive Index of Medium (Water)
        B1 = n1^2;
        B3 = n3^2;
        hbar = 1.054 * 10^-34; % Js, Reduced Planck Constant
        e = -1.602 * 10^-19; % C, Electron Charge
        NA = 6.022 * 10^23; % Avogadro’s number
        E = 78.8; % Relative Permittivity of silica, solvent?
        E0 = 8.85 * 10^-12; % F/m, Permittivity under Vacuum
        c = 3 * 10^8; % m/s, speed of light
        omega = 3.3 * 10^15; % Characteristic Frequency
        
      %% Inputs preprocessing (classic DLVO theory)
        T = TdegC + 273.16; % Convert to absolute temperature
        rho = -0.0025*TdegC^2-0.1249*TdegC+1005.2; % kgm-3 @20MPa and averaged temperature
        mu = 5*10^-9*TdegC^2-3*10^-6*TdegC+0.0005; % Pa s dynamic viscosity
        nu = mu / rho; % m2 s kinetic viscosity
        cion = IS * 1000; % mol m-3, Electrolyte Concentration
        a1nm = d1nm / 2; % Radius of particle 1
        a2nm = d2nm / 2; % Radius of particle 2
        a1 = a1nm * 10^-9; % m, convert to SI unit
        a2 = a2nm * 10^-9; % m, convert to SI unit
        zeta1 = computeZeta(pH, IS)/1000;
        zeta2 = computeZeta(pH, IS)/1000;
        kappa = ((2 * 1^2 * e^2 * cion * NA) / (E * E0 * kB * T))^0.5; % Debye-Hückel parameter for silica
        lambda = c / (pi^2 * omega) * (2 / (B3 * (B1 + B3)))^0.5; % m, Characteristic Wavelength of the Retardation Effect for silica attachment
        
      %% Inputs preprocessing (Soft particle model)
        Hsmin = Hsminnm * 10^-9; % Converting to SI unit
        hcnm(ctd0) = d0nm(ctd0); % Compressing starting distance
        hc(ctd0) = hcnm(ctd0) * 10^-9;
        rhofix = sigma0 / d0(ctd0);
        N0 = rhofix / 1 / e; % Density of free groups
        ydon1 = asinh(1*N0/(2*1*cion*NA));   % Scaled donnan potential
        ydon2 = asinh(1*N0/(1*cion*NA));
        kappa1 = kappa * (cosh(ydon1))^0.5;
        kappa2 = kappa * (cosh(ydon2))^0.5;
        y0 = ydon1 - tanh(ydon1/2); % Scaled unperturbed potential
        phidon1 =  kB * T / (1 * e) * asinh(rhofix/(2*1*e*cion*NA));   % The Donnan potential
        phidon2 =  kB * T / (1 * e) * asinh(rhofix/(2*1*e*cion*NA));
        phio1 = phidon1 - kB * T / (1 * e) * tanh(1*e*phidon1 / (2*kB*T)); % The unperturbed surface potential
        phio2 = phidon2 - kB * T / (1 * e) * tanh(1*e*phidon2 / (2*kB*T));
        phieff1 = 4 * tanh(1*e*phio1/(4*kB*T)); % The effective surface potential
        phieff2 = 4 * tanh(1*e*phio2/(4*kB*T));
        
      %% Separation distance (classic DLVO theory)
        NH = 1000; % Number of nodes
        for j=1:NH % Meshing Hmax = 352081.17, dHmin = 1.39E-15
            dH(j) = 0.0001 / exp((NH + 1 - j)*0.05);
            if j == 1
                H(j) = dH(j);
            else
                H(j) = H(j-1) + dH(j);
            end
            u(j) = H(j) / a1;
        end
        
        %% Separation distance (Soft particle model theory): note - core surface to core surface?
        NHs = NH; % Number of nodes
        for k=1:NHs % Meshing Hmax = 352081.17, dHmin = 1.39E-15
            dHs(k) = 0.0001 / exp((NHs + 1 - k)*0.25);
            if k == 1
                Hs(k) = dHs(k) + Hsmin;
            else
                Hs(k) = Hs(k-1) + dHs(k);
            end
            us(k) = Hs(k) / a1;
        end
        
        %% Solver
        for i=1:NH
            % Attractive potential (Classic)
            A(i) = 0.75 * kB * T * (1 + 2 * kappa * H(i)) * exp(-2 * kappa * H(i)) +...
                3 * hbar * omega / (16 * 2^0.5) * (B1 - B3)^2 / (B1 + B3)^(3/2) *...
                (1 + (H(i) / lambda)^1.185)^(-1/1.185); % Hamaker constants
            VA(i) = -A(i)/6*(2*a1*a2/(H(i)^2 + 2*H(i)*(a1+a2)) + 2*a1*a2/...
                (H(i)^2 + 2*H(i)*(a1+a2) + 4*a1*a2) + log((H(i)^2 + 2*H(i)*...
                (a1+a2))/(H(i)^2 + 2*H(i)*(a1+a2)+4*a1*a2)));
            % Electrostatic potential (classic)
            VE(i) = 2 * pi * E0 * E * 2 * a1 * a2 / (a1 + a2) * zeta1 * zeta2 * log(1 + exp(-kappa*H(i)));
            
            % Attractive potential (Soft particle model)
            As(i) = 0.75 * kB * T * (1 + 2 * kappa * Hs(i)) * exp(-2 * kappa * Hs(i)) +...
                3 * hbar * omega / (16 * 2^0.5) * (B1 - B3)^2 / (B1 + B3)^(3/2) *...
                (1 + (Hs(i) / lambda)^1.185)^(-1/1.185); % Hamaker constants
            VAs(i) = -As(i)/6*(2*a1*a2/(Hs(i)^2 + 2*Hs(i)*(a1+a2)) + 2*a1*a2/...
                (Hs(i)^2 + 2*Hs(i)*(a1+a2) + 4*a1*a2) + log((Hs(i)^2 + 2*Hs(i)*(a1+a2))/...
                (Hs(i)^2 + 2*Hs(i)*(a1+a2)+4*a1*a2)));
            
            % Electrostatic potential (soft particle model)
            if Hs(i) >= 2*d0(ctd0) % Case 0: Approaching
                VEs(i) = 4 * pi * (a1 + d0(ctd0)) * (a2 + d0(ctd0)) / (a1 + d0(ctd0) + a2 + d0(ctd0)) *...
                    E0 * E * phieff1 * phieff2 * exp(-kappa * Hs(i));
                
            else if Hsmin < 2*d0(ctd0) && Hsmin >= d0(ctd0) && Hs(i) < 2*d0(ctd0) && Hs(i) >= Hsmin % Case 1: Interdigitation but no compression -- aggregation happens at H = d0
                    for m = i:NHs
                        G(m) = sinh(kappa1*(Hs(m)-d0(ctd0))) * cosh(kappa2*(Hs(m)/2-d0(ctd0))) -...
                            kappa2/kappa1 * cosh(kappa1*(Hs(m)-d0(ctd0))) * sinh(kappa2*(Hs(m)/2-...
                            d0(ctd0)));
                        if Hs(m)>2*d0(ctd0)
                            G(m)=0;
                            Fpl(m) = 0;
                        else
                            Fpl(m) = 4 * cion * NA * kB * T * (sinh(ydon2/2 - (ydon2-ydon1)/(2*G(m))*...
                                sinh(kappa1*(Hs(m)-d0(ctd0)))*cosh(kappa2*(Hs(m)/2-d0(ctd0))))^2 -...
                                0.25 * (kappa2/kappa)^2 * ((ydon2 - ydon1)/G(m))^2 * sinh(kappa1*(Hs(m)-...
                                d0(ctd0)))^2 * sinh(kappa2*(Hs(m)/2-d0(ctd0)))^2);
                        end
                        
                        if m == i
                            Vpl(m) = 0;
                            Vsp1(m) = 0;
                        else
                            dVpl(m) = 0.5 * (Fpl(m-1)+Fpl(m)) * dHs(m);
                            Vpl(m) = Vpl(m-1) + dVpl(m);
                            dVsp(m) = 2*pi*a1*a2/(a1+a2) * 0.5 * (Vpl(m-1)+Vpl(m)) * dHs(m);
                            Vsp1(m) = Vsp1(m-1) + dVsp(m);
                        end
                    end
                    VEs(i) = Vsp1(NHs);
                    
                else if Hsmin < d0(ctd0) && Hs(i) < d0(ctd0) % Case 2: Compression happens
                        for m = i:NHs
                            Fpl(m) = 2 * cion * NA * kB * T * ((1+(1*N0*d0(ctd0)/(1*cion*NA*Hs(m)))^2)^0.5 - 1);
                            if m == i
                                Vpl(m) = 0;
                                Vsp2(m) = 0;
                            else
                                dVpl(m) = 0.5 * (Fpl(m-1)+Fpl(m)) * dHs(m);
                                Vpl(m) = Vpl(m-1) + dVpl(m);
                                dVsp(m) = 2*pi*a1*a2/(a1+a2) * 0.5 * (Vpl(m-1)+Vpl(m)) * dHs(m);
                                Vsp2(m) = Vsp2(m-1) + dVsp(m);
                            end
                        end
                        VEs(i) = Vsp2(NHs);
                        
                    end
                end
            end
            
            % Total interaction energy
            V(i) = VA(i) + VE(i);
            Vs(i) = VAs(i) + VEs(i);
            
            % Hydrodynamic/viscous effects
            beta(i) = (6*u(i)^2 + 13*u(i) + 2) / (6*u(i)^2 + 4*u(i));
            fW(i) = 2 * beta(i) / (u(i)+2)^2 * exp(V(i)/kB/T);
            betas(i) = (6*us(i)^2 + 13*us(i) + 2) / (6*us(i)^2 + 4*us(i));
            fWs(i) = 2*betas(i)/(us(i)+2)^2*exp(Vs(i)/kB/T);
            
            
            if i == 1
                dWs(i) = 0.5 * fWs(i) * dHs(i);
                Ws(i) = dWs(i);
            else
                dWs(i) = 0.5*(fWs(i) + fWs(i-1))*dHs(i);
                Ws(i) = Ws(i-1) + dWs(i);
            end
            
        end
        
        Wsfinal(ctd0) = Ws(NH);
        W(m,n) = Wsfinal(end);
        if d1nm == 0
            W(m,n) = 0;
        end
        
    end
    
end


end











