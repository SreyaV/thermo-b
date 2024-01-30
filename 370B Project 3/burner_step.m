function [finalState, Tpeak, mdot_mix] = burner(gas, air, mdot_air, P1_air, Tmax,Burner_Pressure_Ratio)

hf = enthalpy_mass(gas);
sf = entroly_mass(gas);
ha = enthalpy_mass(air);
sa = entroly_mass(air);
M_air = meanMolecularWeight(air);
Tpeak = 300;
mdot_mix = 0;

% Choose a small step in fuel.
dfuel = (mdot_air/15)/2000;
for mdot_fuel=dfuel:dfuel:mdot_air
    mdot_mix = mdot_air + mdot_fuel;

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
