function exergy = exergy_mass(gas, flow, T0, P0, environment)

    N = nSpecies(gas);
    M = molecularWeights(gas);
    S = speciesNames(gas);
    CP = chemPotentials(environment);

    if flow
        A = enthalpy_mass(gas) - T0*entropy_mass(gas);
    else
        A = intEnergy_mass(gas) + P0 / density(gas) - T0*entropy_mass(gas);
    end

    tempName = speciesName(gas,1);
    if strcmp(tempName{1},'H2O') && nSpecies(gas) == 1
        % do this for water
        fluid2=Water;
        set(fluid2,'T',T0,'P',P0);
        if flow
            G_tm = enthalpy_mass(fluid2) - T0*entropy_mass(fluid2);
        else
            G_tm = intEnergy_mass(fluid2) + P0 / density(fluid2) - T0*entropy_mass(fluid2);
        end

        Xt = A - G_tm;
        Xc = 0;

        exergy = Xt;
    
    
    elseif strcmp(tempName{1},'C2F4H2') && nSpecies(gas) == 1
        % do this for HFC134a
        fluid2=HFC134a;
        set(fluid2,'T',T0,'P',P0);
        if flow
            G_tm = enthalpy_mass(fluid2) - T0*entropy_mass(fluid2);
        else
            G_tm = intEnergy_mass(fluid2) + P0 / density(fluid2) - T0*entropy_mass(fluid2);
        end

        Xt = A - G_tm;
        Xc = 0;

        exergy = Xt;

    elseif strcmp(tempName{1},'CO2') && nSpecies(gas) == 1
        % do this for CO2
        fluid2=CarbonDioxide;
        tempGas = GRI30;
        set(tempGas,'T',temperature(gas),'P',pressure(gas),'X','CO2:1');
        Xc = exergy_mass(tempGas,flow, T0, P0, environment);
        exergy = Xc;
    
    else
        iN2 = speciesIndex(environment, 'N2');
        iO2 = speciesIndex(environment, 'O2');
        iH2O = speciesIndex(environment, 'H2O');
        iAr = speciesIndex(environment, 'AR');
        iCO2 = speciesIndex(environment, 'CO2');

        gas_IntEx = Solution('gri30.yaml','gri30');
        set(gas_IntEx, 'T', 300, 'P', oneatm, 'X', moleFractions(gas));
        if flow
            A = intEnergy_mass(gas) + P0*1/density(gas) - T0*entropy_mass(gas);
        else
            A = enthalpy_mass(gas) - T0*entropy_mass(gas);
        end
        mass_fracs = massFractions(gas_IntEx);
        mol_fractions = moleFractions(gas);
    
        specs = speciesNames(gas_IntEx);
        Xc_total = 0;
        Xc = 0;
    
        for i=1:length(specs)
            carbons = 0;
            hydrogens = 0;
            nitrogens = 0;
            oxygens = 0;
            argons = 0;
            if mass_fracs(i)>0
                spec = char(S(i))
                carbons = carbons + nAtoms(gas_IntEx, spec, 'C');
                hydrogens = hydrogens + nAtoms(gas_IntEx, spec, 'H');
                nitrogens = nitrogens + nAtoms(gas_IntEx, spec, 'N');
                oxygens = oxygens + nAtoms(gas_IntEx, spec, 'O');
                argons = argons + nAtoms(gas_IntEx, spec, 'Ar');
                
                O2s = (2*carbons + hydrogens/2 - oxygens) /2;
                N2s = nitrogens/2;
                H2Os = hydrogens/2;
                CO2s = carbons;
        
        
                Xc = Xc -  ( CP(iH2O)  * H2Os);
                Xc = Xc - ( CP(iCO2) * CO2s );
                Xc = Xc - ( CP(iN2) * N2s );
                Xc = Xc + ( CP(iO2)  * O2s );
                Xc = Xc;%/M(speciesIndex(gas_IntEx, spec));
                Xc*moleFraction(gas, spec)
                Xc_total = Xc*moleFraction(gas, spec) + Xc_total;
                Xc = 0;
            end
    
        end
    
        Xc_total = Xc_total / meanMolecularWeight(gas)
        
        G0 = Xc_total;
        
        exergy = A + G0;
        
    end
end