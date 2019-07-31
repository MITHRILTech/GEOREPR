function [cHSO4, cSO4] = eq1H2SO4(T, unit, cH, ST, I)

T(T<25) = 25;
I(I<0) = 0;
[Nx,Ny] = size(T);
Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;

% Assume: K1 is infinite.
TdegCdata = [25 30 35 40 43 45 50 60 125 150 175 200 225 250 275 300 325 350];
Sdata = [0.5080 0.5125 0.5176 0.5229 0.5262 0.5282 0.5337 0.5449 0.6422 0.6899 0.7451 0.8097 0.8880 0.9848 1.112 1.287 1.69 2.18];
K20data = [1.028E-2 9.42E-3 6.75E-3 7.25E-3 6.09E-3 6.66E-3 5.31E-3 4.32E-3 7.07E-4 2.75E-4 1.25E-4 5.69E-5 2.65E-5 1.50E-5 4.59E-6 1.94E-6 8.74E-7 3.85E-7];
Akdata = [0.94 0.96 0.98 1.01 1.02 1.03 1.07 1.12 1.42 1.51 1.58 1.65 1.71 1.77 1.77 1.73 1.56 1.34];

for i = 1:Nx
    for j = 1:Ny
        K20(i,j) = interp1(TdegCdata,K20data,TdegC(i,j),'spline');
        S(i,j) = interp1(TdegCdata,Sdata,TdegC(i,j),'spline');
        Ak(i,j) =  interp1(TdegCdata,Akdata,TdegC(i,j),'spline');
    end
end


case1 = (TdegC > 300*ones(size(TdegC)));
case2 = (TdegC <= 300*ones(size(TdegC)));
S = 1.1*S.*case1 + S .* case2;

K2 = 10 .^ (log10(K20) + 4*S.*I.^0.5./(1+Ak.*I.^0.5));


cSO4 = K2 .* ST ./ (cH + K2);
cHSO4 = ST  - cSO4;

end



