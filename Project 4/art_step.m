function gas = atr(w2c, o2c)
    global To Po g_tm mu_o
    g_tm=0;
    gas = Solution('gasification_small.yaml');
    N_species = nSpecies(gas);
    iH2O = speciesIndex(gas, 'H2O');
    iCH4 = speciesIndex(gas, 'CH4');
    iO2 = speciesIndex(gas, 'O2');
    iN2 = speciesIndex(gas, 'N2');
    iH2 = speciesIndex(gas, 'H2');
    iCO2 = speciesIndex(gas, 'CO2');
    iCO = speciesIndex(gas, 'CO');
    iAR = speciesIndex(gas, 'AR');
    M = molecularWeights(gas);
    S = speciesNames(gas);
    Tin = 450+273;
    Pin = 1000000; %Pa
    Po = oneatm;
    To = 298;
    % Air Dead State
    xo=zeros(1, N_species);
    xo(iN2)  = 0.780805;
    xo(iO2)  = 0.209450;
    xo(iAR)  = 0.009340;
    xo(iCO2) = 0.000380;
    w = Water();
    setState_Tsat(w,[To 1]);
    x_H2O = pressure(w)/Po;
    Ntotal = 1/(1-x_H2O);
    for i=1:1:N_species
        xo(i) = xo(i)/Ntotal;
    end
    xo(iH2O) = x_H2O;
    set(gas,'T',To,'P',Po,'X',xo);
    mu_o=chemPotentials(gas);
    majorSpecNames = {'H2O' 'CH4' 'O2' 'N2' 'H2' 'CO2' 'CO'};
    xmix = zeros(1, N_species);
    xmix(iCH4)=1;
    xmix(iO2) = xmix(iCH4)*o2c;
    xmix(iH2O) = xmix(iCH4)*w2c;
    xmix(iN2) = xmix(iCH4)*o2c*3.76;
    set(gas, 'T', Tin, 'P', Pin, 'X', xmix);
    equilibrate(gas,'HP');
end