%% Close Everything
clc;clear;close all;
fprintf('\nClosing everything...\n')

%% Declare Global Variables
global T0 P0 xo mu_o g_tm;                      % Formatted according to Chris' LHV/HHV scripts
g_tm=0;
T0=298.15;
P0=oneatm;

%% Declare Local Variables and Assign Values
fprintf('\nDeclaring Local Variables...\n')
TIT=1410;                                       % Turbine Inlet Temperature from Q2
PR=12.2;                                          % Pressure Ratio from Q2
BurnerPR= 0.95;                                 % Values from Project 1
AirComp_eff=0.86;
FuelComp_eff=0.86;
Turbine_eff=0.82;
HSRG_eff = 0.90;
EconomizerPR = 0.92;
BoilerPR = 0.92;
SuperheaterPR = 0.92;
CondenserPR = 0.92;
SteamTurb_eff = 0.75;
CondensatePump = 0.85;
FeedPump_eff = 0.85;
Condensate_P = 6.8e3; %Pa
TurbineExit_Q = 8.88;
PinchP_TDiff = 20; %K
mdot_air = 144; %kg/s
Air=zeros(1, 4);                                % State variables order-> Enthalpy, Entropy, specific exergy, flow exergy
Fuel=zeros(1, 4);                               % Enthalpy, Entropy, Exergy (Not Specific), Flow Exergy (Not Specific)
Mix=zeros(1, 4);                                % Same as above

%% Declare Fuel
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
set(gas, 'Temperature', T0, 'Pressure', P0, 'X', xfuel);
%% Declare Air
air = Solution('gri30.yaml','gri30');
N = nSpecies(air);
iO2 = speciesIndex(air, 'O2');
iN2 = speciesIndex(air, 'N2');
xair = zeros(1, N);
xair(iO2) = 0.21;
xair(iN2) = 0.79;
set(gas, 'Temperature', T0, 'Pressure', P0, 'X', xair);
%% Declare Water
water = Solution('gri30.yaml','gri30');
N = nSpecies(water);
iH2O = speciesIndex(water, 'H2O');
xwater = zeros(1, N);
xwater(iH2O) = 1;
set(gas, 'Temperature', T0, 'Pressure', P0, 'X', xwater);
%% State 1 -> 2
%[Pret,Tret,finalState,pathStates] = compTurb(fluid,P1,T1,P2,eta_p,nsteps)

%Air
outputs_air = compTurb(air, P0, T0, P0*PR, AirComp_eff, 100);

%Gas
outputs_gas = compTurb(gas, P0, T0, P0*PR*2, AirComp_eff, 100);

%% State 2 -> State 3

%Gas
%outputs_gas = nozzle(gas, 0.5);

%% State 3 -> 4

%Combustor/Mixer

outputs_mix = burner(gas, air, mdot_air, P0*PR, TIT , BurnerPR);
%[finalState, Tpeak, mdot_mix]

%% State 4 -> 5

outputs_mix = compTurb(gas, pressure(gas), temperature(gas), P0, Turbine_eff, 100);
%[Pret,Tret,finalState,pathStates]



