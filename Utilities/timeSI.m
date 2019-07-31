function t = timeSI(tinput, unit)

tfs = strcmpi(unit, 's');
tfmin = strcmpi(unit, 'min');
tfh = strcmpi(unit, 'h');
tfd = strcmpi(unit, 'd');
tfa = strcmpi(unit, 'a');

if tfs == 1
    t = tinput;
elseif tfmin == 1
    t = tinput * 60; 
elseif tfh == 1
    t = tinput * 60 * 60;
elseif tfd == 1
    t = tinput * 60 * 60 * 24;
elseif tfa == 1
    t = tinput * 60 * 60 * 24 * 365;
end

end






