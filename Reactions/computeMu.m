% ref: https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html


function mu = computeMu(T, unit)

Tab = convert2SIunit(T, unit);
TdegC = Tab - 273.15;

% TdegCdata = [0.01	10	20	25	30	40	50	60	70	80	90	100	110	120	140	160	180	200	220	240	260	280	300	320	340	360];
% mudata = [0.0017914	0.001306	0.0010016	0.00089	0.0007972	0.0006527	0.0005465	0.000466	0.0004035	0.000354	0.0003142	0.0002816	0.0002546	0.000232	0.0001966	0.0001704	0.0001504	0.0001346	0.0001218	0.0001111	0.0001018	0.0000936 0.0000859 0.0000783 0.0000703 0.0000603];

p1 =     0.03774;
q1 =       19.67;

mu = (p1) ./ (TdegC + q1);


% mu = griddata(TdegCdata,mudata,TdegC, 'cubic');

end