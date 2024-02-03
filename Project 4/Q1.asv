gas = Solution('gasification_small.yaml','gasification_small');

N = nSpecies(gas);
iH2O = speciesIndex(gas, 'H2O');
iCH4 = speciesIndex(gas, 'CH4');
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
M = molecularWeights(gas);

w2c = 1;
dr = 0.1;
o2c = 0:dr:1;

