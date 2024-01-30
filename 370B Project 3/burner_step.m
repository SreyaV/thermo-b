function [finalState, Tpeak, mdot_mix] = burner(gas, air, mdot_air, P1_air, Tmax,Burner_Pressure_Ratio)

hf = enthalpy_mass(gas);
sf = entropy_mass(gas);
ha = enthalpy_mass(air);
sa = entropy_mass(air);
M_air = meanMolecularWeight(air);
Tpeak = 300;
mdot_mix = 0;
temp   = Solution('gri30.yaml','gri30');
N     = nSpecies(temp);
iCH4  = speciesIndex(gas,'CH4');
iC2H6 = speciesIndex(gas,'C2H6');
iC3H8 = speciesIndex(gas,'C3H8');
iCO2  = speciesIndex(gas,'CO2');
iO2   = speciesIndex(gas,'O2');
iN2   = speciesIndex(gas,'N2');
xfuel = zeros(1,N);
xfuel(iCH4)  = 0.907;
xfuel(iC2H6) = 0.036;
xfuel(iC3H8) = 0.019;
xfuel(iN2)   = 0.018;
xfuel(iCO2)  = 0.010;
xfuel(iO2)   = 0.010;

% Choose a small step in fuel.
dfuel = (mdot_air/15)/2000;
for mdot_fuel=dfuel:dfuel:mdot_air
    mdot_mix = mdot_air + mdot_fuel;
    M_fuel = meanMolecularWeight(gas);
    
    % Mix the air and fuel adiabatically.
    Pm = P1_air;
    hm = (mdot_fuel*hf + mdot_air*ha)/mdot_mix;
    xmix = xfuel*mdot_fuel/M_fuel;
    M_air = meanMolecularWeight(air);
    xmix(iO2)  = xmix(iO2) + 0.21*mdot_air/M_air;
    xmix(iN2)  = xmix(iN2) + 0.79*mdot_air/M_air;
    xmix = xmix/sum(xmix);
    set(gas,'H',hm,'P',Pm,'X',xmix);
    Tm = temperature(gas);
    sm = entropy_mass(gas);
    Sm = sm*mdot_mix;

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

    finalState.T = temperature(gas);
    finalState.P = pressure(gas);
    finalState.H = enthalpy_mass(gas);
    finalState.S = entropy_mass(gas);
end
