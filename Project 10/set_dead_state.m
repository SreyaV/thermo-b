function set_dead_state(gas)
% Set the composition and chemical potentials of the dead state.
% C.F. Edwards, 1-22-11

% Globals for communicating with the main script.
% To Po and RHo have to be set in the calling script.
global To Po RHo xo mu_o

Nsp  = nSpecies(gas);
iO2  = speciesIndex(gas,'O2');
iN2  = speciesIndex(gas,'N2');
iH2O = speciesIndex(gas,'H2O');
iCO2 = speciesIndex(gas,'CO2');
iAR  = speciesIndex(gas,'AR');

% Specify the major components of dry air.  Data are from the Universal
% Industrial Gases website:  http://www.uigi.com/air.html
xo(iN2)  = 0.780805;
xo(iO2)  = 0.209450;
xo(iAR)  = 0.009340;
xo(iCO2) = 0.000380;

% Find the total of the specified species.  Take the remainder that is
% unspecified and put it in as Nitrogen since it is the bulk gas.
xo_specified = 0;
for i=1:1:Nsp
    xo_specified = xo_specified + xo(i);
end
xo_unspecified = 1 - xo_specified;
xo(iN2)  = xo(iN2) + xo_unspecified;

% Add in humidity.  Use a RH spec or molefraction (but be sure it is below
% saturation).
water = Solution('liquidvapor.yaml','water');
setState_Tsat(water,[To 1]);
Psat = pressure(water);
x_H2O_max = Psat/Po
x_H2O = RHo*x_H2O_max
Ntotal = 1/(1-x_H2O);
for i=1:1:Nsp
    xo(i) = xo(i)/Ntotal;
end
xo(iH2O) = x_H2O;
if sum(xo) ~= 1
    disp('Mole fractions do not sum to unity in set_dead_state...')
end

% Set the environmental composition and get the chemical potentials for use
% by the exergy functions.  (They access these through the global
% declaration at the top.)
set(gas,'T',To,'P',Po,'X',xo);
mu_o = chemPotentials(gas);
