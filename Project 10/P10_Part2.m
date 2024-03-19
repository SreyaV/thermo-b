% Part 2 of Project 10

% Global Variables
global e N_A R h F k_B
e = 1.60218e-19;    % Coul/unit-charge
N_A = 6.02214e23;   % particles/mole
R   = 8.31447;      % J/mol-K
h   = 6.62607e-34;  % J-s
F   = e*N_A;        % Coul/mole-of-unit-charge and (J/mol)/(eV/particle)
k_B = R/N_A;        % J/particle-K

%% 

% Set geometries and base variables
channel_width = 10e-3;
channel_height = 5e-3;
channel_length = 0.5;

Tcell = 1000 + 273;
Pcell = 3e5;
inlet_velocity = 1; % m/s
lambda = 2;
excess_air = 1.0;

x_H2 = 0.97;
x_H2O = 1-x_H2;
x_O2 = 1 / 4.76;
x_N2 = 3.76/4.76;

accumulated_current = 0;
diff_current = 0;

voltage = 1; %Change to whatever is desired

% Get a Cantera ideal gas object.  
% Type this for a list of properties:  methods solution -full
gas = Solution('GRI30.yaml');
H2  = speciesIndex(gas,'H2');
H2O = speciesIndex(gas,'H2O');
O2  = speciesIndex(gas,'O2');
N2  = speciesIndex(gas,'N2');
Nsp = nSpecies(gas);

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

% Set the composition and pressure at the cathode:
xcathode = zeros(Nsp,1);
xcathode(O2) = x_O2;
xcathode(N2) = x_N2;
% Set the gas object to the cathode state:
set(gas,'Temperature',Tcell,'Pressure',Pcell,'MoleFractions',xcathode);
% Get the chemical potentials:
mucathode = chemPotentials(gas);    % J/kmol
mucathode = mucathode/1000;         % J/mol

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
muEc_eq = (2/4)*muOc_eq -(1/4)*muO2c_eq;

muEa_eq   = mu_eq(1);
muH2a_eq  = mu_eq(2);
muH2Oa_eq = mu_eq(3);
muO_eq    = mu_eq(4);
muO2c_eq  = mu_eq(5);
muEc_eq   = mu_eq(6);

mu_eq = [muEa_eq muH2a_eq muH2Oa_eq muO_eq muO2c_eq muEc_eq];

%% 

% Evaluate first button cell differential element
[i mu xac] = SOFC_Element_V(voltage,x_eq,mu_eq,Tcell,K,L,ioa,ioc)




