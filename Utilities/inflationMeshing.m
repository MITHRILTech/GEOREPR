function [x, dx, z, dz, zrock, dzrock, t, dt] = inflationMeshing(Nx, Nz, Nzrock, Nt, xmax, zmax, zmaxrock, tmax)

z = zeros(Nz,1);
dz = zeros(Nz,1);
% x = zeros(Nx,1);
% dx = zeros(Nx,1);
t = zeros(Nt,1);
dt = zeros(Nt,1);

zmin = 1e-9;
qz = (zmax/zmin)^(1/(Nz-1));
for i = 1:Nz
    if i == 1
        z(i) = zmin;
        dz(i) =  zmin;
    else
        z(i) = zmin * qz^(i-1);
        dz(i) = z(i) - z(i-1);
    end
end
z = z';

% dz = zmax / Nz;
% z = (dz/2):dz:(zmax-dz/2);
% dz = dz * ones(Nz,1);


dzrock = zmaxrock / Nzrock;
zrock = (dzrock/2):dzrock:(zmaxrock-dzrock/2);




% xmin = 1;
% qx = (xmax / xmin) ^ (1 / (Nx-1));
% for k = 1 : Nx
%     if k == 1
%         x(k) = xmin;
%         dx(k) =  xmin;
%     else
%         x(k) = xmin * qx ^ (k-1);
%         dx(k) = x(k) - x(k-1);
%     end
% end

dx = xmax / Nx;
x = (dx/2):dx:(xmax-dx/2);
% dx = dx * ones(Nx, 1);
x = x';


% dt = tmax/Nt;
% t = (dt/2):dt:(tmax-dt/2);
% dt = tmax/Nt * ones(Nt,1);

tmin = tmax / 10^3;
% tmin = 1;
qt = (tmax / tmin) ^ (1 / (Nt-1));
for j = 1 : Nt
    if j == 1
        t(j) = tmin;
        dt(j) =  tmin;
    else
        t(j) = tmin * qt ^ (j-1);
        dt(j) = t(j) - t(j-1);
    end
end

t = t';

end
