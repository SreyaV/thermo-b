function [t_out, x_out]=oxyHC(gas, fuel_data)
% Preamble %
N_species = nSpecies(gas);
iH2O = speciesIndex(gas, 'H2O');
iCH4 = speciesIndex(gas, 'CH4');
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
iH2 = speciesIndex(gas, 'H2');
iCO2 = speciesIndex(gas, 'CO2');
iCO = speciesIndex(gas, 'CO');
iAR = speciesIndex(gas, 'AR');
iC = elementIndex(gas,'C');
iH = elementIndex(gas,'H');
iO=elementIndex(gas, 'O');
M = molecularWeights(gas);
S = speciesNames(gas);
Tin = 450+273;
Pin = 1000000; %Pa
Po = oneatm;
To = 298;
atom_M = atomicMasses(gas);
% Data %
w2c=fuel_data.w2c;
h2c=fuel_data.h2c;
o2c=fuel_data.o2c;
LHV=fuel_data.LHV;
hfg=fuel_data.hfg;
cp=fuel_data.cp;
mass_fuel=h2c*atom_M(iH)+o2c*atom_M(iO)+atom_M(iC);     
i=1;
for o2_c=0:0.1:1
    m_temp=mass_fuel+w2c*M(iH2O)+3.76*o2_c*M(iN2)+o2_c*M(iO2);
    xmix=zeros(1,N_species);
    xmix([iH2O, iN2, iCO, iO2, iH2])=[w2c, 3.76*o2_c, 1, o2_c + (o2c/2) - 0.5, (h2c/2)];
    if(xmix(iO2)<0)
    xmix(iH2O)=xmix(iH2O)+xmix(iO2)*2;
    xmix(iH2)=xmix(iH2)-xmix(iO2)*2;
    xmix(iO2)=0;
    end
    set(gas,'T',To,'P',Po,'X','CO2:1');
    H_CO2 = enthalpy_mole(gas);
    set(gas,'T',To,'P',Po,'X','H2O:1');
    H_H2O=enthalpy_mole(gas); 
    H_OHC=LHV+hfg+H_CO2/mass_fuel+(h2c/2)*H_H2O/mass_fuel+cp*(425);
    set(gas,'T',Tin,'P',Pin,'X','N2:1');
    H_N2 = enthalpy_mass(gas);
    set(gas,'T',Tin,'P',Pin,'X','O2:1');
    H_O2 = enthalpy_mass(gas);
    set(gas,'T',Tin,'P',Pin,'X','H2O:1');
    H_H2O = enthalpy_mass(gas);
    H_OHC_shifted=(H_OHC*mass_fuel/m_temp)+(H_N2*3.76*o2_c*M(iN2)/m_temp)+(H_O2*o2_c*M(iO2)/m_temp)+(H_H2O*w2c*M(iH2O)/m_temp);
    set(gas,'T',To,'P',Po,'X',xmix);
    equilibrate(gas,'TP');
    set(gas,'H',H_OHC_shifted,'P',Pin);
    equilibrate(gas,'HP');
    t_out(i)=temperature(gas)-273.15;
    x_out(:,i)=moleFractions(gas);
    i=i+1;
end
end