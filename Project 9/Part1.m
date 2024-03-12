%% Button Cell Setup

%Givens

temp = 1000+273;
P   = 100000;

cond_YSZ = 15; % S/m
D_H2H2O  = 3.8378e-3; %m²/s
D_O2N2   = 2.9417e-4;

t_e   = 50e-6; 
t_GDL = 5e-3;

N_H2  = 0.97;
N_H2O = 0.03;

J_cathode = 1000; %A/m^2
J_anode = 100*J_cathode;

% Anode Gas
gas_anode   = Solution('gri30.yaml','gri30');
N     = nSpecies(gas_anode);
iH2   = speciesIndex(gas_anode,'H2');
iH2O   = speciesIndex(gas_anode,'H2O');
M     = molecularWeights(gas_anode);
x_anode = zeros(1,N);
x_anode(iH2O) = N_H2;
x_anode(iH2) = N_H2O;
set(gas_anode,'T',temp,'P',P,'X',x_anode);

% Cathode Gas
gas_cathode   = Solution('gri30.yaml','gri30');
N     = nSpecies(gas_cathode);
iN2   = speciesIndex(gas_cathode,'H2');
iH2O   = speciesIndex(gas_cathode,'H2O');
M     = molecularWeights(gas_cathode);
x_cathode = zeros(1,N);
x_cathode(iH2O) = N_H2;
x_cathode(iH2) = N_H2O;
set(gas_cathode,'T',temp,'P',P,'X',x_cathode);

%% Pass 1 : Cell in Equilibrium

mu = chemicalPotential(gas);


%% 













