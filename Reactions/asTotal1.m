function Astotal = asTotal1(cDSinj, cDSw, Astotal0, T, unit, pH, IS, tinj)

[tinjX, tinjY] = size(tinj);

if tinjX > 1 || tinjY == 1
    tinj0 = [0; tinj(1:end-1)];
elseif tinjX == 1 
    tinj0 = [0 tinj(1:end-1)];
end

dt = tinj - tinj0;
[~, Nz] = size(T);
T = T(:,end);
IS = IS(:,end);
pH = pH(:,end);

% Tab = convert2SIunit(T, unit);
% TdegC = Tab - 273.15;

for i = 1 : length(dt)
    
    if i == 1
        cDSi = cDSw;
        Astotal(i) = Astotal0;
        [SDS, Astotal(i)] = sinkDS(cDSinj, cDSi, Astotal0, T(i), unit, pH(i), IS(i), dt(i));
    else
        [SDS, Astotal(i)] = sinkDS(cDSinj, cDSi, Astotal(i-1), T(i), unit, pH(i), IS(i), dt(i));
    end
        
	
    cDSi = cDSi - SDS*dt(i);
    
    
end

% Astotal = Astotal' * ones(1,Nz);

end