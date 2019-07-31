%Can only be used at pH = [4,8], IS = [0.001, 0.5], assuming temperature
%has no effect on this.

function surfaceCharge = computeSurfacecharge(pH, IS)

ISdata = [0.001; 0.0025; 0.005; 0.007; 0.01; 0.025; 0.05; 0.07; 0.1; 0.25; 0.5; 0.7; 1];
pHdata = [4 6 8];
surfaceChargedata =  [1.39 6.53 50.07; 1.42 8.58 51.38; 1.48 11.30 52.60; 1.51 11.52 55.17; 1.57 11.90 58.20; 1.89 17.06 71.60; 2.43 20.80 82.50; 2.85 23.20 88.20; 3.49 26.50 95.20; 6.69 27.20 107.10; 12.02 27.20 117.20; 16.28 27.21 117.20; 22.67 27.22 117.21];

[Nx, Ny] = size(pH);

surfaceCharge = zeros(Nx, Ny);

pH(isnan(pH)) = 0;
pH(isinf(pH)) = 0;
IS(isnan(IS)) = 0;
IS(isinf(IS)) = 0;


% F = scatteredInterpolant(pHdata,ISdata,surfaceChargedata);

for i = 1:Nx
    for j = 1:Ny
        
        if pH(i,j) > max(pHdata) || IS(i,j) > max(ISdata)
            surfaceCharge(i,j) = griddata(pHdata,ISdata,surfaceChargedata,pH(i,j),IS(i,j), 'nearest');
        else
            surfaceCharge(i,j) = griddata(pHdata,ISdata,surfaceChargedata,pH(i,j),IS(i,j), 'natural');
        end

    end
end




end