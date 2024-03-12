%Button Cell

%Givens

temp = 1000+273;
P = 100000;
cond_YSZ = 15; % S/m
D_H2H2O = 3.8378e-3;
D_O2N2 = 2.9417e-4;

t_e = 50e-6; 
t_GDL = 5e-3;

N_H2 = 0.97;
N_H2O = 0.03;

J_cathode = 1000; %A/m^2
J_anode = 100*J_cathode;

gas_anode   = Solution('gri30.yaml','gri30');

N     = nSpecies(gas);
iH2   = speciesIndex(gas,'H2');
iH2O   = speciesIndex(gas,'H2O');
M     = molecularWeights(gas);
x_anode = zeros(1,N);
x_anode(iH2O) = N_H2;
x_anode(iH2) = N_H2O;
set(gas_anode,'T',temp,'P',P,'X',x_anode);

%% Anode

A_a = 0;

mu = chem_potentials(gas);
mu_H2O = mu(iH2O);
mu_H2 = mu(iH2);
mu_e = 0;

mu_O_YSZ = mu_H2O + mu_e - mu_H2;

%% 
