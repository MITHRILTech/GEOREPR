function Rmd = computeRmd2(cDSi, T, unit, pH, cNaCl)

cDSi(isnan(cDSi)) = 0;
pH(isnan(pH)) = 7;
cNaCl(isnan(cNaCl)) = 0;
h = 0.45;
Tab = convert2SIunit(T, unit);
fNaCl = 1;
aNaCl = cNaCl.*fNaCl; % convert to activity
aNaCl(aNaCl==0) = 0.069;
pHnom = pH + log10(aNaCl./0.069);
St = 10.^(0.0977+75.84./Tab);

% f, f', p137
if pH<=5.97
    f0 = 10.^(pH-7.6);
    fpH = f0./(1+6.2*f0);
    fpHprime = fpH/0.118913;
else
    fpH = 10.^(pH - 7.6 - 2.113 * log10(1+10.^((pH-7.6)/2.113)) - (pH-7.6)./(9.6538 +...
        1.7901 * (pH-7.6) + 4.1811 * (pH-7.6).^2));
    fpHprime = fpH/0.118913;
end

if pHnom<=5.97
    f0nom = 10.^(pHnom-7.6);
    fpHnom = (f0nom./(1+6.2*f0nom));
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

cDS = cDSi;
SSI = cDS./Ce;

if SSI > St
    ff = St.^5 + 5 * St.^4 * (SSI - St);
else
    ff = SSI.^5;
end

Rmd = F .* kOH .* ff .* (1-1./SSI) *0.001/0.01/0.01/60;
Rmd(isnan(Rmd)) = 0;
Rmd(Rmd < 0) = 0;

end