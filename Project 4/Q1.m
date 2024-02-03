gas = Solution('gasification_small.yaml','gasification_small');
mix = Solution('gasification_small.yaml','gasification_small');

N = nSpecies(gas);
iH2O = speciesIndex(gas, 'H2O');
iCH4 = speciesIndex(gas, 'CH4');
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
iH2 = speciesIndex(gas, 'H2');
iCO2 = speciesIndex(gas, 'CO2');
iCO = speciesIndex(gas, 'CO');
M = molecularWeights(gas);
S = speciesNames(gas);

w2c = 1;
dr = 0.01;
o2c = 0:dr:1;

majorSpecNames = {'H2O' 'CH4' 'O2' 'N2' 'H2' 'CO2' 'CO'};
majorSpecies = zeros(length(S), length(o2c));
temperatures = zeros(1, length(o2c));
coldGas_eff = zeros(1, length(o2c));
exergy_eff = zeros(1, length(o2c));
syngasYield = zeros(1, length(o2c));

xgas = zeros(1, N);
xgas(iCH4) = 1;
xgas(iH2O) = w2c;
xgas = xgas./sum(xgas);

Tin = 450+273;
Pin = 1000000; %Pa
Po = oneatm;
To = 298;

set(gas, 'Temperature', Tin, 'Pressure', Pin, 'X', xgas);

xmix = zeros(1, N);

for i=1:length(o2c)
    xmix = xgas;
    xmix(iO2) = xgas(iCH4)*o2c(i);
    xmix(iN2) = xgas(iCH4)*o2c(i)*3.76;
    set(mix, 'Temperature', Tin, 'Pressure', Pin, 'X', xmix);
    
    LHV_initial = LHV_mass(mix);
    exergy_initial = exergy_mass(mix, Po, To);
    
    equilibrate(mix, 'HP');
    
    xmix_normalize = moleFractions(mix);
    for k=1:length(xmix_normalize)
        if xmix_normalize(k)<0.001
            xmix_normalize(k) = 0;
        end
    end
    xmix = xmix_normalize ./ sum(xmix_normalize);
    LHV_final = LHV_mass(mix);
    exergy_final = exergy_mass(mix, Po, To);

    %Major Species
    moleFracs = moleFractions(mix);
    for j=1:length(majorSpecNames)
        index = speciesIndex(mix, majorSpecNames(j));
        if index==0
            majorSpecies(j, i) = 0;
        else
            majorSpecies(j, i) = moleFracs(index);
        end
    end

    %Temperature
    temperatures(i) = temperature(mix);

    %Syngas Yield
    H2Yield = moleFracs(speciesIndex(mix, 'H2'));
    COYield = moleFracs(speciesIndex(mix, 'CO'));
    syngasFrac = 0;
    if H2Yield/COYield < 0.6/0.4
        syngasFrac = H2Yield;
    else
        syngasFrac = COYield;
    end
    syngasYield(i) = syngasFrac/ (w2c/2);

    %Cold-Gas Efficiency
    coldGas_eff(i) = LHV_final/LHV_initial;

    %Exergy Efficiency
    exergy_eff(i) = exergy_initial/exergy_final;

end

figure(1);
hold on;
for k=1:length(majorSpecNames)
    plot(o2c, majorSpecies(k,: ))
end
legend(majorSpecNames)
xlabel("Oxygen/Carbon Molar Feed Ratio (/)")
ylabel("Mole Fraction in Products")
title("Major Species")
hold off;

figure(2);
hold on;
plot(o2c, temperatures)
xlabel("Oxygen/Carbon Molar Feed Ratio (/)")
ylabel("Temperature of Products (K)")
title("Temperature")
hold off;

figure(3);
hold on;
plot(o2c, coldGas_eff)
xlabel("Oxygen/Carbon Molar Feed Ratio (/)")
ylabel("Output/Input LHV Ratio (/)")
title("Cold-Gas Efficiency")
hold off;

figure(4);
hold on;
plot(o2c, exergy_eff)
xlabel("Oxygen/Carbon Molar Feed Ratio (/)")
ylabel("Output/Input Exergy Ratio (/)")
title("Exergy Efficiency")
hold off;

figure(5);
hold on;
plot(o2c, syngasYield)
xlabel("Oxygen/Carbon Molar Feed Ratio (/)")
ylabel("Moles of Syngas Produced per Mole of Input Methane (/)")
title("Syngas Yield")
hold off;



plotfixer_v2()