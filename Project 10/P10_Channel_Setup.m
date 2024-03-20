% Channel setup of Project 10

% Global Variables
global e N_A R h F k_B
e = 1.60218e-19;    % Coul/unit-charge
N_A = 6.02214e23;   % particles/mole
R   = 8.31447;      % J/mol-K
h   = 6.62607e-34;  % J-s
F   = e*N_A;        % Coul/mole-of-unit-charge and (J/mol)/(eV/particle)
k_B = R/N_A;        % J/particle-K

%% VALUES TO SET

% Values you can set!

voltage = 0.5; %Change to whatever is desired
steps = 500; %Change to whatever is desired


%% 

% Set geometries and base variables
channel_width = 10e-3;
channel_height = 5e-3;
channel_length = 0.5;

Tcell = 1000 + 273;
Pcell = 3e5;
inlet_velocity = 0.5; % m/s
lambda = 2;
excess_air = 1.0;

x_H2 = 0.97;
x_H2O = 1-x_H2;
x_O2 = 1 / 4.76;
x_N2 = 3.76/4.76;

accumulated_current = 0;
diff_current = 0;

% Get a Cantera ideal gas object.  
% Type this for a list of properties:  methods solution -full
gas = Solution('GRI30.yaml');
H2  = speciesIndex(gas,'H2');
H2O = speciesIndex(gas,'H2O');
O2  = speciesIndex(gas,'O2');
N2  = speciesIndex(gas,'N2');
Nsp = nSpecies(gas);

% Set the path lengths.
L_YSZ  = 50e-6;                 % m
%L_YSZ  = 0.5e-3;   
L_GDLa = 5e-3;
L_GDLc = 5e-3;
L_Nia  = 0;                     % In case you want the metal too.
L_Nic  = 0;
% Put these in an array for convenient passing.
L = [L_Nia L_GDLa L_GDLa L_YSZ L_GDLc L_Nic];

% Species designation in function Binary_Diffusivity_ijTP:
iN2 = 1; iO2 = 2; iAr = 3; iCO2 = 4; iH2O = 5; iCO = 6; iH2 = 7; iCH4 = 8;

% Set the oxygen ion specific conductivity of the YSZ.
kappa_YSZ = 15;                 % S/m or (A/m2)/(V/m)
Activation_Energy_A=F;
T0=1273;
kappa_YSZ=kappa_YSZ*exp(-(Activation_Energy_A/R)*((1/Tcell)-(1/T0)));

% Since two charges move with each O= ion, the species conductivity
% (product of concentration and molar diffusivity) expressed
% in terms of chemical potential as the driving force is this
% over (zF)^2.
kappa_O = kappa_YSZ/(-2*F)^2 ;   % (mol-O/m2/s)/(J/mol-O/m)
% Set the electron specific conductivity of the nickel.
kappa_Ni = 2e6;                 % S/m or (A/m2)/(V/m)
kappa_E = kappa_Ni/(-1*F)^2 ;    % (mol-E/m2/s)/(J/mol-E/m)
concentration  = Pcell/(R*Tcell);                           % mol/m3
D_O2_N2  = Binary_Diffusivity_ijTP(iO2,iN2,Tcell,Pcell) ;    % m2/s
D_H2_H2O = Binary_Diffusivity_ijTP(iH2,iH2O,Tcell,Pcell);

K(1) = kappa_E;
K(2) = kappa_O;
K(3) = concentration;
K(4) = D_O2_N2;
K(5) = D_H2_H2O;

% Exchange current density of the cathode.
ioco   = 1000;                      % C/s/m2
T_ioco = Tcell;               % K
Ea_ioc = 100e3;                     % J/mol
ioc = ioco*exp(-(Ea_ioc/R)*(1/Tcell - 1/T_ioco));
Anode_Cathode_Ratio = 100;
ioa = Anode_Cathode_Ratio*ioc;

distance_along_channel = zeros(1, steps);
current_density_array = zeros(1, steps);
heat_flux_array = zeros(1, steps);
electrical_power_density = zeros(1, steps);
hydrogen_mole_fractions = zeros(1, steps);
water_mole_fractions = zeros(1, steps);
oxygen_mole_fractions = zeros(1, steps);
equ_electric_potential = zeros(1, steps);
actual_electric_potential = ones(1, steps).*voltage;

dlength = channel_length / steps;
darea = dlength * channel_width;

%% 

% Set arrays to pass to functions

% Set the composition and pressure at the anode:
xanode = zeros(Nsp,1);
xanode(H2O) = x_H2O;
xanode(H2)  = x_H2;
% Set the gas object to the anode state:
set(gas,'Temperature',Tcell,'Pressure',Pcell,'MoleFractions',xanode);
% Get the chemical potentials:
muanode = chemPotentials(gas);      % J/kmol
muanode = muanode/1000;             % J/mol
anode_density = density(gas);
m_anode = meanMolecularWeight(gas);
enthalpy_anode_in = enthalpy_mole(gas);


% Set the composition and pressure at the cathode:
xcathode = zeros(Nsp,1);
xcathode(O2) = x_O2;
xcathode(N2) = x_N2;
% Set the gas object to the cathode state:
set(gas,'Temperature',Tcell,'Pressure',Pcell,'MoleFractions',xcathode);
% Get the chemical potentials:
mucathode = chemPotentials(gas);    % J/kmol
mucathode = mucathode/1000;         % J/mol
air_density = density(gas);
m_air = meanMolecularWeight(gas);
enthalpy_cathode_in = enthalpy_mole(gas);

% Save the needed mole fractions and chemical potentials in their own 
% variable names.
xH2a_eq   = xanode(H2);
xH2Oa_eq  = xanode(H2O);
xO2c_eq   = xcathode(O2);

x_eq = [xH2a_eq xH2Oa_eq xO2c_eq];

muH2a_eq  = muanode(H2);
muH2Oa_eq = muanode(H2O);
muO2c_eq  = mucathode(O2);
% Use the electrochemical potential of an electron in the anode as a
% reference value.  Set this to zero.
muEa_eq = 0;
% The oxygen ion electrochemical potential is determined by the 
% affinity being zero for the anode reaction: H2 + O= -> H2O + 2e-
muO_eq = muH2Oa_eq + 2*muEa_eq - muH2a_eq;
% The oxygen ion has the same electrochemical potential at the cathode when
% in equilibrium.
% The electrochemical potential of an electron at the cathode is given by
% zero affinity for the cathode reaction: O2 + 4e- -> 2O=
muEc_eq = (2/4)*muO_eq -(1/4)*muO2c_eq;

mu_eq = [muEa_eq muH2a_eq muH2Oa_eq muO_eq muO2c_eq muEc_eq];

%ANODE SIDE, NOT AIR
%Getting molar flow rate from anode inlet velocity
%m_air = 28.97; % Molar mass of air, g/mol
A = channel_width*channel_height; % Cross-sectional area of the inlet in m^2
Q = inlet_velocity * A;                            % Volume flow rate in m^3/s
n = anode_density * Q / (0.001*m_anode);     % Moles per second
%kg/m^3 * m^3/s / kg/mol
%Using molar flow rate to find molar flow rate of each gas species
molar_flow_rate_H2 = n*x_H2;
molar_flow_rate_H2O = n*x_H2O;
molar_flow_rate_O2 = molar_flow_rate_H2;
molar_flow_rate_N2 = molar_flow_rate_O2 * (x_N2/x_O2);
initial_H2 = molar_flow_rate_H2;

%TO GET AIR, FIRST FIND HYDROGEN FLOW RATE. KNOW THE STOCHIOMETRY FOR
%HYDROGEN AND OXYGEN. STOCHIOMETRIC WOULD BE THAT OXYGEN IS HALF OF H2, BUT
%SINCE LAMBDA IS 2, THE FLOW RATE OF OXYGEN IN THE CATHODE IS EQUAL TO
%HYDROGEN IN THE ANODE
%BUT IT'S AIR, SO USE THE FLOW RATE OF OXYGEN TO CALCULATE FLOW RATE OF
%NITROGEN TOO, AND ADD.

enthalpy_start_anode = enthalpy_anode_in * (molar_flow_rate_H2 + molar_flow_rate_H2O);
enthalpy_start_cathode = enthalpy_cathode_in * (molar_flow_rate_O2 + molar_flow_rate_N2);


%% 




%% 
% Evaluate first button cell differential element
% walk it to the current
phi_eq = SOFC_Element_icTKL(0,x_eq,mu_eq,Tcell,K,L,ioa,ioc);

walkPotVec = linspace(phi_eq,voltage,100);
i=0;
for walkiter = 1:1:length(walkPotVec)
    [i mu xac delta phi] = SOFC_Element_VTKL(walkPotVec(walkiter),x_eq,mu_eq,Tcell,K,L,ioa,ioc, i);
end

diff_current = i*darea;
accumulated_current = accumulated_current + diff_current;

distance_along_channel(1) = 0;
current_density_array(1) = 0;
heat_flux_array(1) = 0;
electrical_power_density(1) = 0;
hydrogen_mole_fractions(1) = x_H2;
water_mole_fractions(1) = x_H2O;
oxygen_mole_fractions(1) = x_O2;
equ_electric_potential(1) = phi_eq;

%% 

%Evaluate next n-1 elements

iterator = 2;
while iterator <= steps
    
    enthalpy_0 = enthalpy_anode_in* (molar_flow_rate_H2 + molar_flow_rate_H2O);
    enthalpy_0 = enthalpy_cathode_in * (molar_flow_rate_O2 + molar_flow_rate_N2);
    enthalpy_0 = enthalpy_0 / 1000 / darea;

    %Calculate # of electrons used in chemical reaction
    electrons = diff_current / e;
    %Calculate the ratios of anode/cathode gases used/created relative
    %to the electrons
    H2_used = electrons / 2;
    H2O_created = H2_used;
    O2_used = electrons / 4;
    %Use Avogradro's number to calculate the moles / second of each gas
    %that was either used or created in this differential button cell
    %element
    moles_H2_used = H2_used / N_A ;
    moles_H2O_created = H2O_created / N_A;
    moles_O2_used = O2_used / N_A;

    %Update the molar flow rates using the prior values calculated
    molar_flow_rate_O2 = molar_flow_rate_O2 - moles_O2_used;
    molar_flow_rate_N2 = molar_flow_rate_N2; %molar flow rate of n2 should be staying the same but total molar flow rate on cathode side should be decreasing
    %mole fraction of n2 increases. 
    molar_flow_rate_H2 = molar_flow_rate_H2 - moles_H2_used;
    molar_flow_rate_H2O = molar_flow_rate_H2O + moles_H2O_created;
   
    x_H2 = (molar_flow_rate_H2) / (molar_flow_rate_H2 + molar_flow_rate_H2O);
    x_H2O = molar_flow_rate_H2O / (molar_flow_rate_H2O + molar_flow_rate_H2);
    x_O2 = molar_flow_rate_O2 / (molar_flow_rate_O2 + molar_flow_rate_N2);
    x_N2 = molar_flow_rate_N2 / (molar_flow_rate_N2 + molar_flow_rate_O2);
    
    
    %Get the new anode/cathod characterization arrays
    [x_step mu_step enthalpy_anode_out enthalpy_cathode_out] = anode_cathode(x_H2, x_H2O, x_O2, x_N2, Tcell, Pcell);
    %Use the depleted gas arrays to calculate the current density
    %produced by the next differential button cell element
    [i mu xac delta phi] = SOFC_Element_VTKL(voltage,x_step,mu_step,Tcell,K,L,ioa,ioc,i);
    equ_electric_potential(iterator) = SOFC_Element_icTKL(0,x_step,mu_step,Tcell,K,L,ioa,ioc);
    %Calculate the current as current density * area
    diff_current = i*darea;
    accumulated_current = accumulated_current + diff_current;

    enthalpy_1 = enthalpy_anode_out* (molar_flow_rate_H2 + molar_flow_rate_H2O);
    enthalpy_1 = enthalpy_cathode_out * (molar_flow_rate_O2 + molar_flow_rate_N2);
    enthalpy_1 = enthalpy_1 / 1000 / darea;

    distance_along_channel(iterator) = dlength*(iterator);
    current_density_array(iterator) = i;
    electrical_power_density(iterator) = i*voltage;
    heat_flux_array(iterator) = enthalpy_0 - enthalpy_1 - i*voltage;
    hydrogen_mole_fractions(iterator) = x_H2;
    water_mole_fractions(iterator) = x_H2O;
    oxygen_mole_fractions(iterator) = x_O2;
    %equ_electric_potential(iterator) = phi;

    enthalpy_anode_in = enthalpy_anode_out;
    enthalpy_cathode_in = enthalpy_cathode_out;


    iterator = iterator + 1;
    %disp(i);
end

%For Graph 1
accumulated_current
power = accumulated_current*voltage;
fuel_utilization = 1 - molar_flow_rate_H2 / initial_H2;
enthalpy_out_anode = enthalpy_anode_out * (molar_flow_rate_H2 + molar_flow_rate_H2O);
enthalpy_out_cathode = enthalpy_cathode_out * (molar_flow_rate_O2 + molar_flow_rate_N2);
heat_transfer = (enthalpy_start_anode+enthalpy_start_cathode)/1000 - (enthalpy_out_anode+enthalpy_out_cathode)/1000 - power;

%enthalpy of reactants going in, H2 and O2, N2, H2O
%enthalpy of products coming out, same things coming out
%enthalpy of reactants coming in = enthalpy coming out + output power +
%heat transfer
%use enthalpy function from cantera, its by kmol, divide by 1000 to get
%per mole, convert to mass


%For Graph 2
%plot all the following:

distance_along_channel;
current_density_array;
heat_flux_array ;
electrical_power_density ;
hydrogen_mole_fractions ;
water_mole_fractions;
oxygen_mole_fractions;
equ_electric_potential;
actual_electric_potential;

%% 
close all;
figure(1)
hold on;
plot(distance_along_channel, current_density_array./1000)
plot(distance_along_channel, heat_flux_array.*-1.9/1000)
plot(distance_along_channel, electrical_power_density./1000)
xlabel("Distance Along Channel (m)")
title("1000 C, 3 bar, 0.5 V")
legend("Current Density (kA/m^2)", "Heat Flux (kW/m^2)", "Elec. Power Density (kW/m^2)")
hold off;

plotfixer()

%% 

figure(2)
hold on;
plot(distance_along_channel, hydrogen_mole_fractions)
plot(distance_along_channel, oxygen_mole_fractions)
plot(distance_along_channel, water_mole_fractions)
plot(distance_along_channel, equ_electric_potential)
plot(distance_along_channel, actual_electric_potential, 'k--')
xlabel("Distance Along Channel (m)")
title("1000 C, 3 bar, 0.61 V")
legend("xH2 (anode)", "xO2 (cathode)", "xH2O (anode)", "\Delta \phi equil. (V)", "\Delta \phi cell (V)")
hold off;

plotfixer()

