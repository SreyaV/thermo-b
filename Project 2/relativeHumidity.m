%Code

rh = [linspace(0.0001, 0.1, 17) linspace(0.1, 1, 27)]
satP = 0.0313 *oneatm;
xi_water = satP/oneatm .* rh;


exergy_co = zeros(1, length(xi_water));
exergy_ch4 = zeros(1, length(xi_water));
exergy_h2 = zeros(1, length(xi_water));

%% 

species = 'CH4';

gas = Solution('gri30.yaml','gri30');
N = nSpecies(gas);

iElem = speciesIndex(gas, species);
xgas = zeros(1, N);
xgas(iElem) = 1;
set(gas, 'Temperature', 300, 'Pressure', oneatm, 'X', xgas);

for j=1:length(xi_water)
    
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
    xenv(iH2O) = xi_water(j);
    xenv(iAr) = 0.009015;
    xenv(iCO2) = 0.000397;
    xenv = xenv./sum(xenv);
    set(environment, 'T', 300, 'P', oneatm, 'X', xenv);
    CP = chemPotentials(environment);
    mol_fractions = moleFractions(environment);

%%%%
    gas_IntEx = Solution('gri30.yaml','gri30');
    set(gas_IntEx, 'T', 300, 'P', oneatm, 'X', moleFractions(gas));
    Gtm = intEnergy_mass(gas) + oneatm*1/density(gas) - 300*entropy_mass(gas)

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
    
    exergy_ch4(j) = Gtm + Xc_total;

end

%% 

figure(1);
hold on;
plot(rh.*100, exergy_co./exergy_co(end))
plot(rh.*100, exergy_ch4./exergy_ch4(end))
plot(rh.*100, exergy_h2./exergy_h2(end))
xlabel("Relative Humidity of Environment (%)")
ylabel("Chemical Exergy Relative to 100% Humidity")
legend('CO', 'CH4', 'H2')
scale = axis;
axis([0 100 0.98 1.1])
hold off;
plotfixer_v2()
%% 

%Function
