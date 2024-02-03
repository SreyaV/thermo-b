% Calculation of the effect of varying dead state humidity on the chemical
% exergy of carbon monoxide, methane, and hydrogen.
% C.F. Edwards, 1/8/11

clear all
format compact
fprintf('\n***************************************************************\n')

% Make global variables to communicate the dead state.
global To Po mu_o

% Make a gas object to use for analysis.
gas       = Solution('gri30.yaml','gri30');
Nspecies  = nSpecies(gas);
Nelements = nElements(gas);

% Define indices to find the environmental species in the array.
iCO2 = speciesIndex(gas,'CO2');
iH2O = speciesIndex(gas,'H2O');
iN2  = speciesIndex(gas,'N2');
iO2  = speciesIndex(gas,'O2');
iAR  = speciesIndex(gas,'AR');

% Find the fuels too.
iCH4 = speciesIndex(gas,'CH4');
iCO  = speciesIndex(gas,'CO');
iH2 = speciesIndex(gas,'H2');

% Set the state of the environment.  Fix a model for dry air and then add
% the humidity in afterwards.
To = 25+273.15
Po = 101325

% Specify the major components of air but without the humidity.
xo = zeros(1,Nspecies);
xo(iN2)  = 0.757223;
xo(iO2)  = 0.202157;
xo(iAR)  = 0.009015;
xo(iCO2) = 0.000397;
% Renormalize to make dry air molefractions.
xo = xo/sum(xo);
xN2dry  = xo(iN2)
xO2dry  = xo(iO2)
xARdry  = xo(iAR)
xCO2dry = xo(iCO2)

% Save a copy of this composition.  We will need to reset to this in our
% humidity loops.
xoo = xo;

% Add in humidity.  Find the max water content.
water = Solution('liquidvapor.yaml','water');
Psat = satPressure(water,To)
x_H2O_max = Psat/Po;

% Run a loop through the list of fuels to compare.
iFuel = [iCH4 iCO iH2];
for k=1:1:length(iFuel)
    k
    % Run through a loop of RH.
    j = 1;
    for RH_o=1:-.0001:.0001
        % Put the humidity into the dry air.
        x_H2O = RH_o*x_H2O_max;
        Ntotal = 1/(1-x_H2O);
        
        % Don't forget to normalize.
        xo = xoo;
        for i=1:1:Nspecies
            xo(i) = xo(i)/Ntotal;
        end
        xo(iH2O) = x_H2O;
        xo = xo/sum(xo);
        
        % Set the environmental composition and get the chemical
        % potentials.
        set(gas,'T',To,'P',Po,'X',xo);
        mu_o = chemPotentials(gas);

        % Set the sample gas compositions.
        x = zeros(1,Nspecies);
        x(iFuel(k)) = 1;
        set(gas,'T',To,'P',Po,'X',x);

        % Get and save the chemical exergy.
        RH(j)    = RH_o;
        X_C(k,j) = exergy_mass(gas);
        
        j = j+1;
    end
end

figure(1)
plot(100*RH,X_C(2,:)/X_C(2,1),'r')
hold on
plot(100*RH,X_C(1,:)/X_C(1,1),'g')
plot(100*RH,X_C(3,:)/X_C(3,1),'b')
legend('CO','CH_4','H_2')
xlabel('Relative Humidty of Environment (%)')
ylabel('Exergy at RH/Saturation Exergy') 
hold off
axis([-5 100 .99 1.1])
plotfixer

% Make a plot of mass-specific exergy so you can check the calculations.

figure(2)
plot(100*RH,X_C(2,:)/1e6,'r')
hold on
plot(100*RH,X_C(1,:)/1e6,'g')
plot(100*RH,X_C(3,:)/1e6,'b')
legend('CO','CH_4','H_2')
xlabel('Relative Humidty of Environment (%)')
ylabel('Chemical Exergy (MJ/kg)') 
hold off
scale = axis;
axis([-5 100 scale(3) scale(4)])
plotfixer
