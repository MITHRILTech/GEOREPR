function Ew = relativepermittivityWater(T, unit)

% at 10 MPa
[Nx,Ny] = size(T);
Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;

TdegCdata = [0 25 50 75 100 125 150 175 200 225 250 275 300];
Ewdata = [88 79 70 62.5 55.5  49.5 44.5 39.5 35 31 27.5 24 20.5];

Ew = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny
        Ew(i,j) =  interp1(TdegCdata,Ewdata,TdegC(i,j),'spline');
    end
end

end