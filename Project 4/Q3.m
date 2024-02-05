%% Q3 : Generate H2 Yield and Equilibrium Species arrays

clc;close all;clear;
global To Po g_tm mu_o

g_tm=0;

gas = Solution('gasification_small.yaml');
gas_no_methane = Solution('gasification_no_methane.yaml');

%Different for gas and gas without methane
N_species = nSpecies(gas);
iH2O = speciesIndex(gas, 'H2O');
iCH4 = speciesIndex(gas, 'CH4');
iO2 = speciesIndex(gas, 'O2');
iN2 = speciesIndex(gas, 'N2');
iN2m = speciesIndex(gas_no_methane, 'N2');
iH2 = speciesIndex(gas, 'H2');
iNH3 = speciesIndex(gas, 'NH3');
iCO2 = speciesIndex(gas, 'CO2');
iCO2m = speciesIndex(gas_no_methane, 'CO2');
iCO = speciesIndex(gas, 'CO');
iCOm = speciesIndex(gas_no_methane, 'CO');
iAR = speciesIndex(gas, 'AR');
M = molecularWeights(gas);
S = speciesNames(gas);

Tin = 450+273;
Pin = 10e5; %Pa
Po = oneatm;
To = 298;

% Air Dead State : needed just for exergy (not here)
xo=zeros(1, N_species);
xo(iN2)  = 0.780805;
xo(iO2)  = 0.209450;
xo(iAR)  = 0.009340;
xo(iCO2) = 0.000380;
set(gas,'T',To,'P',Po,'X',xo);

mu_o  =  chemPotentials(gas);
majorSpecNames = {'H2O' 'CH4' 'N2' 'H2' 'CO2' 'CO'};

res = 0.01;
w2c = res:res:3;
o2c = res/10:res:1;

% Looping through H2O/C ratio and O2/C ratio
for j = 1:length(w2c)
   
    for i = 1:length(o2c) 
        

        %%% ATR

        %Initial gas composition : CH4, H2O, O2, N2
        xmix = zeros(1, N_species);
        xmix(iCH4) =1;
        xmix(iO2)  = xmix(iCH4)*o2c(i);
        xmix(iH2O) = xmix(iCH4)*w2c(j);
        xmix(iN2)  = xmix(iCH4)*o2c(i)*3.76;
        xmix = xmix/sum(xmix);
        initial_CH4 = xmix(iCH4);
        set(gas, 'T', Tin, 'P', Pin, 'X', xmix);

        % Oxidize 
        equilibrate(gas,'HP');
        xmix_postATR = moleFractions(gas);

        % Number of moles using C conservation, initial number of moles:1
        moles_number = xmix(iCH4)/ (xmix_postATR(iCO2) +xmix_postATR(iCO) + xmix_postATR(iCH4)) ;

        %%% Shift Reactor
        
        % First store and then remove NH3 and CH4
        xshift = moleFractions(gas);
        nCH4 = xshift(iCH4)*moles_number;
        nNH3 = xshift(iNH3)*moles_number;
        moles_number_in = moles_number - nCH4 -nNH3;
        xshift(iNH3) = []; % Remove NH3 first because it is last in the list
        xshift(iCH4) = [];
        xshift = xshift/sum(xshift);
        set(gas_no_methane, 'T', Tin-200, 'P', pressure(gas), 'X', xshift);

        % Water shift reaction
        equilibrate(gas_no_methane, 'TP');
        xmix_postReactor = moleFractions(gas_no_methane);
        
        % Compute number of moles in outlet 
        moles_number_out =  (xshift(iCO2m) +xshift(iCOm))/ (xmix_postReactor(iCO2m) +xmix_postReactor(iCOm)) *moles_number_in ;
        xshift = moleFractions(gas_no_methane)*moles_number_out;

        % Add stored NH3 and CH4
        xmix_final = [xshift(1:iCH4-1); nCH4; xshift(iCH4:iNH3-1); nNH3; xshift(iNH3:end)];
        set(gas,'T', temperature(gas_no_methane),'P',pressure(gas_no_methane),'X',xmix_final);
        
        % Final quantities
        moleFracs = moleFractions(gas);
        moles_out = moles_number_out + nNH3 +nCH4;
        H2Yield(j, i) = moleFracs(iH2)*moles_out/initial_CH4; %j=w2c and i=o2c 
        

        %Major Species @w2c = 1
        if w2c(j) ==1
            for l=1:length(majorSpecNames) %l=name; i=o2c
                index = speciesIndex(gas, majorSpecNames(l));
                if index==0
                    majorSpecies(l, i) = 0;
                else
                    majorSpecies(l, i) = moleFracs(index);
                end
            end
        end

    end
end

%% Plotting

%Shifted Molar H2 Yield
clf
figure(1)
[x, y]=contour(o2c,w2c,H2Yield);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('Shifted Molar H2 Yield (H2/CH4)')
plotfixer
legend('off')

% Shifted Equilibrium Mole Fraction
majorSpecNames = {'H2O' 'CH4' 'N2' 'H2' 'CO2' 'CO'};
figure(2)
hold on;
for k=1:length(majorSpecNames)
    plot(o2c, majorSpecies(k,: ))
end
legend(majorSpecNames)
xlabel("Oxygen/Carbon Molar Feed Ratio")
ylabel("Shifted Equilibrium Mole Fraction (H2O/C : 1.0)")
hold off;
plotfixer