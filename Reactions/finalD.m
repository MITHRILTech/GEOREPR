% T = [298 303 353 358 368 371 373 433 473 473 523 543 568 568 613 613 623];
% d = [8 7 5 6 7 7 8 15 10 12 15 15 36 64 88 105 150];

function dfinal = finalD(T, unit)

Tab = convert2SIunit(T, unit);
% TdegC = Tab - 273.15;

x = Tab;

a =   1.706e-05;
b =       25.42;
c =       7.707;

dfinal =  a*exp(0.001*b*x)+c;
       
% General model:
%      f(x) = a*exp(0.001*b*x)+c
% Coefficients (with 95% confidence bounds):
%        a =   1.706e-05  (-7.368e-05, 0.0001078)
%        b =       25.42  (16.83, 34.01)
%        c =       7.707  (0.6549, 14.76)
% 
% Goodness of fit:
%   SSE: 1339
%   R-square: 0.9548
%   Adjusted R-square: 0.9483
%   RMSE: 9.78