function [interState,finalState, Tpeak, mdot_mix,mdot_fuel,Tm] = burner(gas, air, mdot_air, P1_air, Tmax,Burner_Pressure_Ratio)

iCH4 = speciesIndex(gas, 'CH4');
iC2H6 = speciesIndex(gas, 'C2H6');
iC3H8 = speciesIndex(gas, 'C3H8');
iCO2 = speciesIndex(gas, 'CO2');
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');

hf = enthalpy_mass(gas);
sf = entropy_mass(gas);
ha = enthalpy_mass(air);
sa = entropy_mass(air);
xair = moleFractions(air);
M_air = meanMolecularWeight(air);
Tpeak = 300;
mdot_mix = 0;
xfuel = moleFractions(gas);
M_fuel = meanMolecularWeight(gas);

% Choose a small step in fuel.
dfuel = (mdot_air/15)/2000;
for mdot_fuel=dfuel:dfuel:mdot_air
    mdot_mix = mdot_air + mdot_fuel;

    % Mix the air and fuel adiabatically.
    Pm = P1_air;
    hm = (mdot_fuel*hf + mdot_air*ha)/mdot_mix;
    xmix = xfuel*mdot_fuel/M_fuel+xair*mdot_air/M_air;
    xmix = xmix/sum(xmix);
    set(gas,'H',hm,'P',Pm,'X',xmix);
    Tm = temperature(gas);
    sm = entropy_mass(gas);
    Sm = sm*mdot_mix;
    interState.T = temperature(gas);
    interState.P = pressure(gas);
    interState.h = enthalpy_mass(gas);
    interState.s = entropy_mass(gas);
    interState.x = exergy_mass(gas);
    interState.xf = flowExergy_mass(gas);


    % Burn the mixture adiabatically.
    Pm_burn = Pm*Burner_Pressure_Ratio;
    set(gas,'H',hm,'P',Pm_burn);
    equilibrate(gas,'HP');
    Tm_burn = temperature(gas);
    
    % Look for a simple crossover since the steps are small and fast.
    if(Tm_burn >= Tmax)
        Tpeak = Tm_burn
        break
    end
end
mdot_fuel
finalState.T = temperature(gas);
finalState.P = pressure(gas);
finalState.h = enthalpy_mass(gas);
finalState.s = entropy_mass(gas);
finalState.x = exergy_mass(gas);
finalState.xf = flowExergy_mass(gas);

end
