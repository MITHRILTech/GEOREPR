% function [cHSO4, cSO4] = eqH2SO4(T, unit, cT, cH, I)
% 
% Tab = convert2SIunit(T, unit);
% TdegC = Tab - 273.15;
% 
% % Assume: K1 is infinite.
% TdegCdata = [25 30 35 40 43 45 50 60 125 150 175 200 225 250 275 300 325 350];
% Sdata = [0.5080 0.5125 0.5176 0.5229 0.5262 0.5282 0.5337 0.5449 0.6422 0.6899 0.7451 0.8097 0.8880 0.9848 1.112 1.287 1.69 2.18];
% K20data = [1.028E-2 9.42E-3 6.75E-3 7.25E-3 6.09E-3 6.66E-3 5.31E-3 4.32E-3 7.07E-4 2.75E-4 1.25E-4 5.69E-5 2.65E-5 1.50E-5 4.59E-6 1.94E-6 8.74E-7 3.85E-7];
% Akdata = [0.94 0.96 0.98 1.01 1.02 1.03 1.07 1.12 1.42 1.51 1.58 1.65 1.71 1.77 1.77 1.73 1.56 1.34];
% K20 = interp1(TdegCdata,K20data,TdegC,'spline');
% S = interp1(TdegCdata,Sdata,TdegC,'spline');
% Ak =  interp1(TdegCdata,Akdata,TdegC,'spline');
% 
% if TdegC > 300
%     S = 1.1 * S;
% end
% K2 = 10 ^ (log10(K20) + 4*S*I^0.5/(1+Ak*I^0.5));
% 
% % a = 1;
% % b = (cH2SO40 + cH0 + cSO40 + K2);
% % c = (cSO40 * (cH2SO40 + cH0) - (cH2SO40 + cHSO40) * K2);
% % 
% % if ((-b + (b^2 - 4*a*c)^0.5) / (2 * a)) >= 0 
% %     x = (-b + (b^2 - 4*a*c)^0.5) / (2 * a);
% % else
% %     x = (-b - (b^2 - 4*a*c)^0.5) / (2 * a);
% % end
% % 
% % cH2SO4 = 0;
% % cHSO4 = cH2SO40 + cHSO40 - x;
% % cSO4 = cSO40 + x;
% % cH = cH2SO40 + cH0 + x;
% 
% cHSO4 = cH * cT / (K2 + cH);
% cSO4 = cT - cHSO4;
% 
% end

function [cH2SO4, cH, cHSO4, cSO4] = eqH2SO4(T, unit, cH0, cH2SO40, cHSO40, cSO40, I)

Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;

% Assume: K1 is infinite.
TdegCdata = [25 30 35 40 43 45 50 60 125 150 175 200 225 250 275 300 325 350];
Sdata = [0.5080 0.5125 0.5176 0.5229 0.5262 0.5282 0.5337 0.5449 0.6422 0.6899 0.7451 0.8097 0.8880 0.9848 1.112 1.287 1.69 2.18];
K20data = [1.028E-2 9.42E-3 6.75E-3 7.25E-3 6.09E-3 6.66E-3 5.31E-3 4.32E-3 7.07E-4 2.75E-4 1.25E-4 5.69E-5 2.65E-5 1.50E-5 4.59E-6 1.94E-6 8.74E-7 3.85E-7];
Akdata = [0.94 0.96 0.98 1.01 1.02 1.03 1.07 1.12 1.42 1.51 1.58 1.65 1.71 1.77 1.77 1.73 1.56 1.34];
K20 = interp1(TdegCdata,K20data,TdegC,'spline');
S = interp1(TdegCdata,Sdata,TdegC,'spline');
Ak =  interp1(TdegCdata,Akdata,TdegC,'spline');
if TdegC > 300
    S = 1.1 * S;
end
K2 = 10 ^ (log10(K20) + 4*S*I^0.5/(1+Ak*I^0.5));

a = 1;
b = (cH2SO40 + cH0 + cSO40 + K2);
c = (cSO40 * (cH2SO40 + cH0) - (cH2SO40 + cHSO40) * K2);

x(1) = (-b + (b^2 - 4*a*c)^0.5) / (2 * a);
x(2) = (-b - (b^2 - 4*a*c)^0.5) / (2 * a);
x(x < 0) = [];


cH2SO4 = 0;
cHSO4 = cH2SO40 + cHSO40 - min(x);
cSO4 = cSO40 + min(x);
cH = cH2SO40 + cH0 + min(x);

end





    
    




