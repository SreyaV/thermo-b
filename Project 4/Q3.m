clc;close all;clear;
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
res=0.01;
j=1;
for w2c=0:res:3
    i=1;
    for o2c=0:res*10:1 
        xmix = zeros(1, N_species);
        xmix(iCH4)=1;
        xmix(iO2) = xmix(iCH4)*o2c;
        xmix(iH2O) = xmix(iCH4)*w2c;
        xmix(iN2) = xmix(iCH4)*o2c*3.76;
        set(gas, 'T', Tin, 'P', Pin, 'X', xmix);
        LHV_initial = LHV_mass(gas);
        massFracs_in = massFractions(gas);
        exergy_initial = flowExergy_mass(gas);
        equilibrate(gas,'HP');

        %%%End ATR
        %%%Start Shift Reactor
    
        gas_shift = Solution('gasification_no_methane.yaml');
        moleFracs = moleFractions(gas);
        COFrac = moleFracs(iCO);
        H2OFrac = 1-COFrac;
        xgas_shift = zeros(1, N_species);
        xgas_shift(iCO) = COFrac;
        xgas_shift(iH2O) = H2OFrac;
        set(gas_shift, 'T', 250+273, 'P', Pin, 'X', xgas_shift);
        equilibrate(gas_shift, 'TP');
        
        xpost_shift = zeros(1, N_species);
        
        xpost_shift(iCO) = 


        LHV_final = LHV_mass(gas);
        massFracs_out  = massFractions(gas);
        moleFracs = moleFractions(gas);
        exergy_final = flowExergy_mass(gas);
        temperatures(j, i) = temperature(gas)-273.15;
        exergy_eff(j, i) = exergy_final/exergy_initial;
        coldGas_eff(j, i) = LHV_final/LHV_initial;
        H2Yield(j, i) = moleFracs(iH2);
        COYield(j, i) = moleFracs(iCO);
        CO2Yield(j, i) = moleFracs(iCO2);
        O2Yield(j, i) = moleFracs(iO2);
        CH4Yield(j, i) = moleFracs(iCH4);
        N2Yield(j, i) = moleFracs(iN2);
        H2OYield(j, i) = moleFracs(iH2O);
        syngasYield(j, i) =(massFracs_out(iCO)/M(iCO)+massFracs_out(iH2)/M(iH2))/(massFracs_in(iCH4)/M(iCH4));
        i=i+1;
    end
    j=j+1;
end