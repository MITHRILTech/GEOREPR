function x = lengthSI(xinput, unit)

tfkm = strcmpi(unit, 'km');
tfm = strcmpi(unit, 'm');
tfmm = strcmpi(unit, 'mm');
tfum = strcmpi(unit, 'um');
tfnm = strcmpi(unit, 'nm');

if tfm == 1
    x = xinput;
elseif tfkm == 1
    x = xinput * 1000;
elseif tfmm == 1
    x = xinput * 0.001;
elseif tfum == 1
    x = xinput * 0.001 * 0.001;
elseif tfnm == 1
    x = xinput * 0.001 * 0.001 * 0.001;
end
    
end

