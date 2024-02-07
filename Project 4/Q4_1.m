clc;clear;close all;
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
iC = elementIndex(gas,'C');
iH = elementIndex(gas,'H');
iO=elementIndex(gas, 'O');
M = molecularWeights(gas);
S = speciesNames(gas);
Tin = 450+273;
Pin = 1000000; %Pa
Po = oneatm;
To = 298;
atomic_M = atomicMasses(gas);
M_CH4=meanMolecularWeight(Methane);
LHV_CH4=50050*10^3;         %Table A-27
CP_CH4=2.20*10^3;           %Table A-27
w2c=1;
i=1;
for o2c=0:0.01:1
    xmix = zeros(1, N_species);
    xmix(iCH4)=1;
    xmix(iO2) = xmix(iCH4)*o2c;
    xmix(iH2O) = xmix(iCH4)*w2c;
    xmix(iN2) = xmix(iCH4)*o2c*3.76;
    set(gas, 'T', Tin, 'P', Pin, 'X', xmix);
    equilibrate(gas,'HP');
    moleFracs(:,i) = moleFractions(gas);
    temperatures(i) = temperature(gas)-273.15;
    massFracs_indirect=(xmix.*M')./(sum(xmix.*M'));
    set(gas,'T',To,'P',Po,'X','CO2:1');
    H_CO2=enthalpy_mole(gas);
    set(gas,'T',To,'P',Po,'X','H2O:1');
    H_H2O=enthalpy_mole(gas);
    H_CH4=LHV_CH4+CP_CH4*(425)+H_CO2/M_CH4+2*H_H2O/M_CH4;
    set(gas,'T',Tin,'P',Pin,'X','N2:1');
    H_N2=enthalpy_mass(gas);
    set(gas,'T',Tin,'P',Pin,'X','O2:1');
    H_O2=enthalpy_mass(gas);
    set(gas,'T',Tin,'P',Pin,'X','H2O:1');
    H_H2O=enthalpy_mass(gas);
    H_CH4_indirect=H_CH4*massFracs_indirect(iCH4)+H_O2*massFracs_indirect(iO2)+H_N2*massFracs_indirect(iN2)+H_H2O*massFracs_indirect(iH2O);
    xmix_indirect=xmix;
    xmix_indirect(iCH4)=0;
    xmix_indirect(iCO)=xmix(iCH4);
    xmix_indirect(iH2)=xmix(iCH4)*2;
    xmix_indirect(iO2)=xmix(iO2)-xmix(iCH4)/2;
    if(xmix_indirect(iO2)<0)
    xmix_indirect(iH2O)=xmix_indirect(iH2O)+xmix_indirect(iO2)*2;
    xmix_indirect(iH2)=xmix_indirect(iH2)-xmix_indirect(iO2)*2;
    xmix_indirect(iO2)=0;
    end
    set(gas,'T',To,'P',Po,'X',xmix_indirect);
    equilibrate(gas,'TP');
    set(gas,'H',H_CH4_indirect,'P',Pin);
    equilibrate(gas,'HP');
    moleFracs_indirect(:,i)=moleFractions(gas);
    temperatures_indirect(i)=temperature(gas)-273.15;
    i=i+1;
end
o2c=0:0.01:1;
figure(1)
hold on
plot(o2c,temperatures_indirect,'--')
plot(o2c,temperatures, 'r')
yline(450, '--k')
legend('Direct','Indirect')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Equilibrium Temperature')
text(.1,1400,sprintf('Methane Autothermal Reforming'))
text(.1,1300,'Conditions: 10 bar, 450 deg preheat')
text(.1,1200,'Water/Methane Mole Ratio: 1')
hold off
plotfixer
figure(2)
hold on
plot(o2c,moleFracs(iH2,:),'b')
plot(o2c,moleFracs(iN2,:),'y')
plot(o2c,moleFracs(iCO,:),'m')
plot(o2c,moleFracs(iCO2,:),'c')
plot(o2c,moleFracs(iH2O,:),'g')
plot(o2c,moleFracs(iCH4,:),'r')
yline(0,'k')
yline(0,'k--')
plot(o2c,moleFracs_indirect(iH2,:),'b--')
plot(o2c,moleFracs_indirect(iN2,:),'y--')
plot(o2c,moleFracs_indirect(iCO,:),'m--')
plot(o2c,moleFracs_indirect(iCO2,:),'c--')
plot(o2c,moleFracs_indirect(iH2O,:),'g--')
plot(o2c,moleFracs_indirect(iCH4,:),'r--')
legend('H2','N2','CO','CO2','H2O','CH4','Direct','Indirect')
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Equilibrium Mole Fraction')
text(.1,0.45,sprintf('Methane Autothermal Reforming'))
text(.1,0.40,'Conditions: 10 bar, 450 deg preheat')
text(.1,0.35,'Water/Methane Mole Ratio: 1')
hold off
plotfixer