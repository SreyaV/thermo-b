% Make a table of exergy and heating values for ideal gases in GRI30.
% C.F. Edwards, 1/8/11

clear all
format compact
fprintf('\n***********************************************************\n')

% Make the environmental (dead) state available to all functions.
global To Po xo mu_o

% Make a Cantera gas object and get indices for species of interest.
gas  = Solution('gri30.yaml','gri30');
N    = nSpecies(gas);
iO2  = speciesIndex(gas,'O2');
iN2  = speciesIndex(gas,'N2');
iCO2 = speciesIndex(gas,'CO2');
iH2O = speciesIndex(gas,'H2O');
iAR  = speciesIndex(gas,'AR');

iH2   = speciesIndex(gas,'H2');
iCO   = speciesIndex(gas,'CO');
iCH4  = speciesIndex(gas,'CH4');
iC3H8 = speciesIndex(gas,'C3H8');
iC2H6 = speciesIndex(gas,'C2H6');

% Set up a water object to work with for saturation.
water = Solution('liquidvapor.yaml','water');
% Set the environmental state.
To = 25+273.15
Po = 101325
Psat = satPressure(water,To)
xsat = Psat/Po
xo = zeros(1,N);
xo(iN2)  = 0.757223;
xo(iO2)  = 0.202157;
xo(iH2O) = 0.031208;
xo(iAR)  = 0.009015;
xo(iCO2) = 0.000397;
% Check for sum to unity.
Sum_xo = sum(xo)
if Sum_xo ~= 1
    disp('Mole fractions dont sum to unity.')
end

% Get the chemical potentials.  These will be used by the exergy function
% and are accessed through global storage.
set(gas,'T',To,'P',Po,'X',xo);
mu_o = chemPotentials(gas);

% Do some species/state calculations.
xgas = zeros(1,N);
xgas(iH2) = 1;
set(gas,'T',25+273.15,'P',101325,'X',xgas);
LHV_H2 = LHV_mass(gas);
HHV_H2 = HHV_mass(gas);
xc_H2 = exergy_mass(gas);
fc_H2 = flowExergy_mass(gas);

xgas = zeros(1,N);
xgas(iCO) = 1;
set(gas,'T',25+273.15,'P',101325,'X',xgas);
LHV_CO = LHV_mass(gas);
HHV_CO = HHV_mass(gas);
xc_CO = exergy_mass(gas);
fc_CO = flowExergy_mass(gas);

xgas = zeros(1,N);
xgas(iCH4) = 1;
set(gas,'T',25+273.15,'P',101325,'X',xgas);
LHV_CH4 = LHV_mass(gas);
HHV_CH4 = HHV_mass(gas);
xc_CH4 = exergy_mass(gas);
fc_CH4 = flowExergy_mass(gas);

xgas = zeros(1,N);
xgas(iC3H8) = 1;
set(gas,'T',25+273.15,'P',101325,'X',xgas);
LHV_C3H8 = LHV_mass(gas);
HHV_C3H8 = HHV_mass(gas);
xc_C3H8 = exergy_mass(gas);
fc_C3H8 = flowExergy_mass(gas);

xgas = zeros(1,N);
xgas(iN2) = 1;
set(gas,'T',25+273.15,'P',101325,'X',xgas);
LHV_N2 = LHV_mass(gas);
HHV_N2 = HHV_mass(gas);
xc_N2 = exergy_mass(gas);
fc_N2 = flowExergy_mass(gas);

xgas = zeros(1,N);
xgas(iO2) = 1;
set(gas,'T',25+273.15,'P',101325,'X',xgas);
LHV_O2 = LHV_mass(gas);
HHV_O2 = HHV_mass(gas);
xc_O2 = exergy_mass(gas);
fc_O2 = flowExergy_mass(gas);

xgas = zeros(1,N);
xgas(iCO2) = 1;
set(gas,'T',25+273.15,'P',101325,'X',xgas);
LHV_CO2 = LHV_mass(gas);
HHV_CO2 = HHV_mass(gas);
xc_CO2 = exergy_mass(gas);
fc_CO2 = flowExergy_mass(gas);

xgas = zeros(1,N);
xgas(iN2)   = 1;
xgas(iCH4)  = 0.907;
xgas(iC2H6) = 0.036;
xgas(iC3H8) = 0.019;
xgas(iN2)   = 0.018;
xgas(iCO2)  = 0.010;
xgas(iO2)   = 0.010;
set(gas,'T',25+273.15,'P',101325,'X',xgas);
LHV_NatGas = LHV_mass(gas);
HHV_NatGas = HHV_mass(gas);
xc_NatGas = exergy_mass(gas);
fc_NatGas = flowExergy_mass(gas);

xgas = zeros(1,N);
xgas(iCO) = 0.4;
xgas(iH2) = 0.6;
set(gas,'T',25+273.15,'P',10e6,'X',xgas);
LHV_SynGas = LHV_mass(gas);
HHV_SynGas = HHV_mass(gas);
xc_SynGas = exergy_mass(gas);
fc_SynGas = flowExergy_mass(gas);

xgas = zeros(1,N);
xgas(iN2) = 0.79;
xgas(iO2) = 0.21;
set(gas,'T',25+273.15,'P',101325,'X',xgas);
LHV_EngrAir = LHV_mass(gas);
HHV_EngrAir = HHV_mass(gas);
xc_EngrAir = exergy_mass(gas);
fc_EngrAir = flowExergy_mass(gas);
set(gas,'T',25+273.15,'P',10e6,'X',xgas);
LHV_CompAir = LHV_mass(gas);
HHV_CompAir = HHV_mass(gas);
xc_CompAir = exergy_mass(gas);
fc_CompAir = flowExergy_mass(gas);
set(gas,'T',0+273.15,'P',101325,'X',xgas);
LHV_ColdAir = LHV_mass(gas);
HHV_ColdAir = HHV_mass(gas);
xc_ColdAir = exergy_mass(gas);
fc_ColdAir = flowExergy_mass(gas);
set(gas,'T',650+273.15,'P',101325,'X',xgas);
LHV_WarmAir = LHV_mass(gas);
HHV_WarmAir = HHV_mass(gas);
xc_WarmAir = exergy_mass(gas);
fc_WarmAir = flowExergy_mass(gas);

% Make a table.
fprintf('\n')
fprintf('Gas                LHV (J/kg)     HHV (J/kg)   Exergy (J/kg)  Flow Exergy (J/kg)\n')
fprintf('--------------------------------------------------------------------------------\n')
fprintf('Hydrogen           %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_H2,HHV_H2,xc_H2,fc_H2)
fprintf('Carbon Monoxide    %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_CO,HHV_CO,xc_CO,fc_CO)
fprintf('Methane            %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_CH4,HHV_CH4,xc_CH4,fc_CH4)
fprintf('Propane            %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_C3H8,HHV_C3H8,xc_C3H8,fc_C3H8)
fprintf('Nitrogen           %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_N2,HHV_N2,xc_N2,fc_N2)
fprintf('Oxygen             %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_O2,HHV_O2,xc_O2,fc_O2)
fprintf('Carbon Dioxide     %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_CO2,HHV_CO2,xc_CO2,fc_CO2)
fprintf('Natural Gas        %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_NatGas,HHV_NatGas,xc_NatGas,fc_NatGas)
fprintf('Simple Syngas      %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_SynGas,HHV_SynGas,xc_SynGas,fc_SynGas)
fprintf('Engineering Air    %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_EngrAir,HHV_EngrAir,xc_EngrAir,fc_EngrAir)
fprintf('Compressed E-Air   %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_CompAir,HHV_CompAir,xc_CompAir,fc_CompAir)
fprintf('Cold E-Air         %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_ColdAir,HHV_ColdAir,xc_ColdAir,fc_ColdAir)
fprintf('Warm E-Air         %10.3e     %10.3e     %10.3e     %10.3e\n',LHV_WarmAir,HHV_WarmAir,xc_WarmAir,fc_WarmAir)
fprintf('\n')

% Make a table for the TM 39 class.
set(gas,'T',20+273.15,'P',2e6,'X',xgas);
LHV_CompAir = LHV_mass(gas);
HHV_CompAir = HHV_mass(gas);
xc_CompAir = exergy_mass(gas);
fc_CompAir = flowExergy_mass(gas);
fprintf('\n')
fprintf('Gas                T (C)    P (bar)    LHV (MJ/kg)    HHV (MJ/kg)    Exergy (MJ/kg)\n')
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('Hydrogen           %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_H2/1e6,HHV_H2/1e6,xc_H2/1e6)
fprintf('Carbon Monoxide    %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_CO/1e6,HHV_CO/1e6,xc_CO/1e6)
fprintf('Methane            %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_CH4/1e6,HHV_CH4/1e6,xc_CH4/1e6)
fprintf('Propane            %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_C3H8/1e6,HHV_C3H8/1e6,xc_C3H8/1e6)
fprintf('Nitrogen           %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_N2/1e6,HHV_N2/1e6,xc_N2/1e6)
fprintf('Oxygen             %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_O2/1e6,HHV_O2/1e6,xc_O2/1e6)
fprintf('Carbon Dioxide     %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_CO2/1e6,HHV_CO2/1e6,xc_CO2/1e6)
fprintf('Natural Gas        %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_NatGas/1e6,HHV_NatGas/1e6,xc_NatGas/1e6)
fprintf('Simple Syngas      %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_SynGas/1e6,HHV_SynGas/1e6,xc_SynGas/1e6)
fprintf('Engineering Air    %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',To-273.15,Po/1e5,LHV_EngrAir/1e6,HHV_EngrAir/1e6,xc_EngrAir/1e6)
fprintf('Compressed E-Air   %5.1f     %5.1f    %10.3f     %10.3f     %10.3f\n',20,2e6/1e5,LHV_CompAir/1e6,HHV_CompAir/1e6,xc_CompAir/1e6)
fprintf('Cold E-Air         %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',0,Po/1e5,LHV_ColdAir/1e6,HHV_ColdAir/1e6,xc_ColdAir/1e6)
fprintf('Warm E-Air         %5.1f     %5.3f    %10.3f     %10.3f     %10.3f\n',650,Po/1e5,LHV_WarmAir/1e6,HHV_WarmAir/1e6,xc_WarmAir/1e6)
fprintf('\n')
