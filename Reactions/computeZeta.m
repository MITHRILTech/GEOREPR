function zetaPotential = computeZeta(pH, IS)

ISdata = [0.001; 0.0025; 0.005; 0.007; 0.01; 0.025; 0.05; 0.07; 0.1];
pHdata = [2.6 4 6 8];
zetadata =  [-1 -24 -36 -45; 1 -22 -33 -42; 3 -15 -31 -37; 2.5 -12 -24 -33; 2 -9 -22 -26; 2.1 -5 -21 -23; 2.2 -6 -19 -19.5; 3 -6 -17.5 -18; 5 -7 -16 -20];

zetaPotential = griddata(pHdata, ISdata, zetadata, pH, IS, 'cubic');

end