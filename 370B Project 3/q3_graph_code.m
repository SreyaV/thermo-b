%% Close Everything
clc;clear;close all;
fprintf('\nClosing everything...\n')

%% Declare Global Variables
global T0 P0 xo mu_o g_tm;                      % Formatted according to Chris' LHV/HHV scripts

g_tm =0;
T0   =298.15;
P0   =oneatm;

%% Declare Local Variables and Assign Values
fprintf('\nDeclaring Local Variables...\n')
TIT           = 1410;                     % Turbine Inlet Temperature from Q2
PR            = 12.2;                     % Pressure Ratio from Q2
BurnerPR      = 0.95;                     % Values from Project 1
AirComp_eff   = 0.86;
FuelComp_eff  = 0.86;
Turbine_eff   = 0.82;
HSRG_eff      = 0.90;
EconomizerPR  = 0.92;
BoilerPR      = 0.92;
SuperheaterPR = 0.92;
CondenserPR   = 0.92;
SteamTurb_eff = 0.75;
CondPump_eff  = 0.85;
FeedPump_eff  = 0.85;
Condensate_P  = 6.8e3;  %Pa
TurbineExit_Q = 0.88;
PinchP_TDiff  = 20;     %K
mdot_air      = 144;    %kg/s
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
fuelStates.f1.T = T0;
fuelStates.f1.P = P0;
fuelStates.f1.h = enthalpy_mass(gas);
fuelStates.f1.s = entropy_mass(gas);

%% Declare Air
air = Solution('gri30.yaml','gri30');
N = nSpecies(air);
iO2 = speciesIndex(air, 'O2');
iN2 = speciesIndex(air, 'N2');
xair = zeros(1, N);
xair(iO2) = 0.21;
xair(iN2) = 0.79;
set(air, 'Temperature', T0, 'Pressure', P0, 'X', xair);
airStates.a1.T = T0;
airStates.a1.P = P0;
airStates.a1.h = enthalpy_mass(air);
airStates.a1.s = entropy_mass(air);

%% Declare Water
water = Water;
set(water, 'Temperature', T0, 'Pressure', P0);
Tw1 = temperature(water);
hw1 = enthalpy_mass(water);
sw1 = entropy_mass(water);

%% State 1af -> 2af

%Air
[Pret_1a_2a,Tret_1a_2a,airStates.a2,path_1a_2a] = compTurb(air, P0, T0, P0*PR, AirComp_eff, 50);

%Gas
[Pret_1f_2f,Tret_1f_2f,fuelStates.f2,path_1f_2f] = compTurb(gas, P0, T0, P0*PR*2, AirComp_eff, 50);

%% State 2af -> State 3m

%Gas
nozzle(gas, 0.5);
hm3  = enthalpy_mass(gas);

%% State 3m -> 4m

%Combustor/Mixer
[State_4m, Tpeak, mdot_mix,mdot_fuel,Tm3] = burner(gas, air, mdot_air, P0*PR, TIT , BurnerPR);

%% State 4m -> 5m

[Pret_4m_5m,Tret_4m_5m,State_5m,path_4m_5m] = compTurb(gas, pressure(gas), temperature(gas), P0, Turbine_eff, 100);


%% State 1w -> 2w : Find first P water after feed pump

pressure_ratio_fpump = 125; 
P2w = P0*pressure_ratio_fpump; 
[Pret1w_2w,Tret_1w_2w,waterStates.w1,path_1w_2w] = pump(water,P0,P2w,FeedPump_eff,100);



%% State 2w -> 5w && 5m -> 8m : Find first mdot_w before HX

mix = gas; % Define the mix coming out of the gas turbine

mdot_w = 25.5; % kg/s
[Tpinch,waterStates,mixStates] = HRSG(mix,water,mdot_mix,mdot_w,HSRG_eff,EconomizerPR,BoilerPR,SuperheaterPR);

% Careful here because HRSG does not set "water"
set(water,'T',waterStates.w4.T,'P', waterStates.w4.P);
Tpinch

%% State 5w -> 6w

[Pret_5w_6w,Tret_5w_6w,waterStates.w5,pathStates_5w_6w] = compTurb(water,pressure(water),temperature(water),Condensate_P,Turbine_eff,10);
Quality = vaporFraction(water)

%% State 6w -> 7w

waterStates.w6 =condenser(water, Condensate_P, CondenserPR, TurbineExit_Q);


%% State 7w -> 8w

[Pret_7w_8w,Tret_7w_8w,waterStates.w7,pathStates_7w_8w] = pump(water,pressure(water),P0,CondPump_eff,100);

%% States

% Make arrays of plotting points.
    
Tw1 = T0;
Tw2 = waterStates.w1.T;
Tw3 = waterStates.w2.T;
Tw4 = waterStates.w3.T;
Tw5 = waterStates.w4.T;
Tw6 = waterStates.w5.T;
Tw7 = waterStates.w6.T;
Tw8 = waterStates.w7.T;


hw1 = hw1;
hw2 = waterStates.w1.h;
hw3 = waterStates.w2.h;
hw4 = waterStates.w3.h;
hw5 = waterStates.w4.h;
hw6 = waterStates.w5.h;
hw7 = waterStates.w6.h;
hw8 = waterStates.w7.h;


Tf1 = fuelStates.f1.T;
Tf2 = fuelStates.f2.T;


hf1 = fuelStates.f1.h;
hf2 = fuelStates.f2.h;


Ta1 = airStates.a1.T;
Ta2 = airStates.a2.T;


ha1 = airStates.a1.h;
ha2 = airStates.a2.h;


Tm3 = Tm3;
Tm4 = State_4m.T;
Tm5 = mixStates.m1.T;
Tm6 = mixStates.m2.T;
Tm7 = mixStates.m3.T;
Tm8 = mixStates.m4.T;


hm3 = (hf2*mdot_fuel+ha2*mdot_air)/mdot_mix;
hm4 = State_4m.h;
hm5 = mixStates.m1.h;
hm6 = mixStates.m2.h;
hm7 = mixStates.m3.h;
hm8 = mixStates.m4.h;


%% Total work

% Gas turbine
W_gross_gas_turbine = mdot_mix*(hm4 - hm5);
W_fuel_compressor   = mdot_fuel*(hf2 - hf1);
W_air_compressor    = mdot_air*(ha2 - ha1);
W_net_gas_turbine = W_gross_gas_turbine - W_fuel_compressor - W_air_compressor;

%Steam turbine
W_gross_steam_turbine = mdot_w*(hw5 - hw6);
W_feed_pump   = mdot_w*(hw2 - hw1);
W_cond_pump   = mdot_w*(hw8 - hw7);
W_net_steam_turbine = W_gross_steam_turbine - W_feed_pump - W_cond_pump;

%% PLOTS


%% T-h plot
Th_plot = true;

if Th_plot 
    clf
    q3_Th_plot
end
