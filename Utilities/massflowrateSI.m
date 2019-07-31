function mDot = massflowrateSI(mInput, unit, Tinj)

rho = -0.0025*Tinj^2-0.1249*Tinj+1005.2; %kgm-3 @20MPa and averaged temperature

% mass unit: g, kg, ton;
% volume unit: m3, dm3(L), cm3(mL);
% time unit: s, min, h;

tf1 = strcmpi(unit, 'g/s');
tf2 = strcmpi(unit, 'kg/s');
tf3 = strcmpi(unit, 'ton/s');
tf4 = strcmpi(unit, 'g/min');
tf5 = strcmpi(unit, 'kg/min');
tf6 = strcmpi(unit, 'ton/min');
tf7 = strcmpi(unit, 'g/h');
tf8 = strcmpi(unit, 'kg/h');
tf9 = strcmpi(unit, 'ton/h');
tf10 = strcmpi(unit, 'm3/s');
tf11 = strcmpi(unit, 'dm3/s');
tf12 = strcmpi(unit, 'cm3/s');
tf13 = strcmpi(unit, 'm3/min');
tf14 = strcmpi(unit, 'dm3/min');
tf15 = strcmpi(unit, 'cm3/min');
tf16 = strcmpi(unit, 'm3/h');
tf17 = strcmpi(unit, 'dm3/h');
tf18 = strcmpi(unit, 'cm3/h');


if tf1 == 1 % g/s
    mDot = mInput / 1000;
elseif tf2 == 1 % kg/s
    mDot = mInput;
elseif tf3 == 1 % ton/s
    mDot = mInput * 1000;
elseif tf4 == 1 % g/min
    mDot = mInput / 1000 / 60;
elseif tf5 == 1 % kg/min
    mDot = mInput / 60;
elseif tf6 == 1 % ton/min
    mDot = mInput * 1000 / 60;
elseif tf7 == 1 % g/h
    mDot = mInput / 1000 / 3600;
elseif tf8 == 1 % kg/h
    mDot = mInput / 3600;
elseif tf9 == 1 % ton/h
    mDot = mInput * 1000 / 3600;
elseif tf10 == 1 % m3/s
    mDot =  mInput * rho; 
elseif tf11 == 1 % dm3/s
    mDot = mInput /1000 * rho; 
elseif tf12 == 1 % cm3/s
    mDot = mInput /1000/1000 * rho; 
elseif tf13 == 1 % m3/min
    mDot = mInput * rho / 60; 
elseif tf14 == 1 % dm3/min
    mDot = mInput /1000 * rho / 60; 
elseif tf15 == 1 % cm3/min
    mDot = mInput /1000/1000 * rho / 60; 
elseif tf16 == 1 % m3/h
    mDot = mInput * rho / 3600; 
elseif tf17 == 1 % dm3/h
    mDot = mInput /1000 * rho / 3600; 
elseif tf18 == 1 % cm3/h
    mDot = mInput /1000/1000 * rho / 3600; 
end

% if tfmkg == 0
%     mDot = mInput * rho / 3600; % t/h to kg/s
% else
%     mDot = mInput; % kg/s
% end

end
    