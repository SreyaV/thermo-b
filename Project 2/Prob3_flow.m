%% Code

species = ["H2", "CO", "CH4", "C3H8", "N2", "O2", "CO2", "Natural Gas", "Simplified Syngas", "Engineering Air", "Compressed Engineering Air", "Cold Engineering Air", "Warm Engineering Air"];

int_ex_flow_values = zeros(1, length(species));


% gas = Solution('gri30.yaml','gri30');
% N = nSpecies(gas);
% 
% iElem = speciesIndex(gas, species);
% xgas = zeros(1, N);
% xgas(iElem) = 1;
% set(gas, 'Temperature', 300, 'Pressure', oneatm, 'X', xgas);

% int_ex = exergy_mass(gas)

%Tracker
tracker = 1;

%Single Species
single_species = ["H2", "CO", "CH4", "C3H8", "N2", "O2", "CO2"];
single_species_cell = cellstr(single_species);

for i=1:length(single_species)
    gas = Solution('gri30.yaml','gri30');
    N = nSpecies(gas);
    iElem = speciesIndex(gas, single_species_cell(i));
    xfuel = zeros(1, N);
    xfuel(iElem) = 1;
    single_species_cell(i)
    set(gas, 'Temperature', 298, 'Pressure', oneatm, 'X', xfuel);

    int_ex = exergy_mass(gas)
    int_ex_flow_values(tracker) = int_ex;

    tracker = tracker + 1;
end


%Natural Gas

gas = Solution('gri30.yaml','gri30');
N = nSpecies(gas);
iCH4 = speciesIndex(gas, 'CH4');
iC2H6 = speciesIndex(gas, 'C2H6');
iC3H8 = speciesIndex(gas, 'C3H8');
iCO2 = speciesIndex(gas, 'CO2');
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
M = molecularWeights(gas);

xfuel = zeros(1, N);
xfuel(iCH4) = 0.907;
xfuel(iC2H6) = 0.036;
xfuel(iC3H8) = 0.019;
xfuel(iCO2) = 0.010;
xfuel(iN2) = 0.018;
xfuel(iO2) = 0.010;
set(gas, 'Temperature', 298, 'Pressure', oneatm, 'X', xfuel);

int_ex = exergy_mass(gas)
int_ex_flow_values(tracker) = int_ex;
tracker = tracker + 1;


% Simplified Syngas
gas = Solution('gri30.yaml','gri30');
N = nSpecies(gas);
iCO = speciesIndex(gas, 'CO');
iH2 = speciesIndex(gas, 'H2');
xfuel = zeros(1, N);
xfuel(iCO) = 0.4;
xfuel(iH2) = 0.6;
set(gas, 'Temperature', 298, 'Pressure', 10e6, 'X', xfuel);

int_ex = exergy_mass(gas)
int_ex_flow_values(tracker) = int_ex;
tracker = tracker + 1;

%Engineering Air
gas = Solution('gri30.yaml','gri30');
N = nSpecies(gas);
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
xfuel = zeros(1, N);
xfuel(iO2) = 0.21;
xfuel(iN2) = 0.79;
set(gas, 'Temperature', 298, 'Pressure', oneatm, 'X', xfuel);

int_ex = exergy_mass(gas)
int_ex_flow_values(tracker) = int_ex;
tracker = tracker + 1;

%Compressed Engineering Air
gas = Solution('gri30.yaml','gri30');
N = nSpecies(gas);
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
xfuel = zeros(1, N);
xfuel(iO2) = 0.21;
xfuel(iN2) = 0.79;
set(gas, 'Temperature', 298, 'Pressure', 10e6, 'X', xfuel);

int_ex = exergy_mass(gas)
int_ex_flow_values(tracker) = int_ex;
tracker = tracker + 1;

%Cold Engineering Air
gas = Solution('gri30.yaml','gri30');
N = nSpecies(gas);
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
xfuel = zeros(1, N);
xfuel(iO2) = 0.21;
xfuel(iN2) = 0.79;
set(gas, 'Temperature', 273, 'Pressure', oneatm, 'X', xfuel);

int_ex = exergy_mass(gas)
int_ex_flow_values(tracker) = int_ex;
tracker = tracker + 1;

%Warm Engineering Air
gas = Solution('gri30.yaml','gri30');
N = nSpecies(gas);
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
xfuel = zeros(1, N);
xfuel(iO2) = 0.21;
xfuel(iN2) = 0.79;
set(gas, 'Temperature', 923, 'Pressure', oneatm, 'X', xfuel);

int_ex = exergy_mass(gas)
int_ex_flow_values(tracker) = int_ex;
tracker = tracker + 1;

%% Function

function internal_exergy = exergy_mass(gas)
    
    N = nSpecies(gas);
    M = molecularWeights(gas);
    S = speciesNames(gas);

    %set up environment
    environment = Solution('gri30.yaml','gri30');
    N = nSpecies(environment);
    iN2 = speciesIndex(environment, 'N2');
    iO2 = speciesIndex(environment, 'O2');
    iH2O = speciesIndex(environment, 'H2O');
    iAr = speciesIndex(environment, 'AR');
    iCO2 = speciesIndex(environment, 'CO2');
    xenv = zeros(1, N);
    xenv(iN2) = 0.757223;
    xenv(iO2) = 0.202157;
    xenv(iH2O) = 0.031208;
    xenv(iAr) = 0.009015;
    xenv(iCO2) = 0.000397;
    set(environment, 'T', 298, 'P', oneatm, 'X', xenv);
    CP = chemPotentials(environment);
    mol_fractions = moleFractions(environment)

    %%%%%%%%%%%%%%%%%%%%%%%%%

    gas_IntEx = Solution('gri30.yaml','gri30');
    set(gas_IntEx, 'T', 298, 'P', oneatm, 'X', moleFractions(gas));
    Gtm = enthalpy_mass(gas) - 298*entropy_mass(gas)

    mass_fracs = massFractions(gas_IntEx);
    mol_fractions = moleFractions(gas);

    specs = speciesNames(gas_IntEx)
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
    Gtm


    % Total exergy is the sum of thermo-mechanical and chemical exergy
    internal_exergy = Gtm + Xc_total;
end
