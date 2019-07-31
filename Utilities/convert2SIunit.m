function x = convert2SIunit(input, unit)

% Comfirmed, can be used for matrix

lfkm = strcmpi(unit, 'km');
lfm = strcmpi(unit, 'm');
lfmm = strcmpi(unit, 'mm');
lfum = strcmpi(unit, 'um');
lfnm = strcmpi(unit, 'nm');
tfs = strcmpi(unit, 's');
tfmin = strcmpi(unit, 'min');
tfh = strcmpi(unit, 'h');
tfd = strcmpi(unit, 'd');
tfa = strcmpi(unit, 'a');
TfdegC = strcmpi(unit, 'degC');
TfK = strcmpi(unit, 'K');


if lfm == 1
    x = input;
elseif lfkm == 1
    x = input * 1000;
elseif lfmm == 1
    x = input * 0.001;
elseif lfum == 1
    x = input * 0.001 * 0.001;
elseif lfnm == 1
    x = input * 0.001 * 0.001 * 0.001;
elseif tfs == 1
    x = input;
elseif tfmin == 1
    x = input * 60; 
elseif tfh == 1
    x = input * 60 * 60;
elseif tfd == 1
    x = input * 60 * 60 * 24;
elseif tfa == 1
    x = input * 60 * 60 * 24 * 365;
elseif TfdegC == 1
    x = input + 273.15;
elseif TfK == 1
    x = input;
end
    
end
