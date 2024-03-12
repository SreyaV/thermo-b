%% Button Cell Setup
clear, close all

%Givens
disp("Setting up..")
F = 96485.33212; % Faraday cst (Coulomb /mole)
R = 8.3144; % Perfect gas cst (SI units)
temp = 1000+273; % K
P   = 1e5; % Pa

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
M_gas_anode     = molecularWeights(gas_anode);
x_anode = zeros(1,N);
x_anode(iH2) = N_H2;
x_anode(iH2O) = N_H2O;
set(gas_anode,'T',temp,'P',P,'X',x_anode);

% Cathode Gas
gas_cathode   = Solution('gri30.yaml','gri30');
N     = nSpecies(gas_cathode);
iN2   = speciesIndex(gas_cathode,'N2');
iO2   = speciesIndex(gas_cathode,'O2');
x_cathode = zeros(1,N);
x_cathode(iO2) = 0.21;
x_cathode(iN2) = 0.79;
set(gas_cathode,'T',temp,'P',P,'X',x_cathode);
M_gas_cathode     = molecularWeights(gas_cathode);

%% Pass 1 : Cell in Equilibrium

disp("Equilibrium cell... Electrochemical potential is (J/mol):")
mu_anode_eq = chemPotentials(gas_anode)/1000; % J/mol
mu_cathode_eq = chemPotentials(gas_cathode)/1000;

mu_electron_anode_eq = 0; %Convention
mu_O2ion_YSZ_anode_eq  = mu_anode_eq(iH2O) - mu_anode_eq(iH2) + 2*mu_electron_anode_eq;
mu_O2ion_YSZ_cathode_eq = mu_O2ion_YSZ_anode_eq;

mu_electron_cathode_eq   = 0.5 * mu_O2ion_YSZ_cathode_eq - 0.25*mu_cathode_eq(iO2);

delta_phi_eq = 0.5/F * ( mu_anode_eq(iH2) + 0.5 * mu_cathode_eq(iH2O) -mu_anode_eq(iH2O))


%% Pass 2

disp("Pass 2: Out of Equilibrium ...")
dj = 0.1;
jmax = 50;
current_densities = 1:dj:jmax;

for i=1  %:length(current_densities)
    j = 1%current_densities(i)*1000; % Current density A/m²
    v = (1/2/F)*j; % Reaction velocity per area : mol/s/m²

    % Fluxes per area: anode
    J_H2_anode  = v;
    J_H2O_anode = v;
    J_electron  = 2*v;

    % Chemical potential differences (PER AREA) @anode
    c = P/ (R*temp); % Total concentration of species
    delta_mu_e_anode = J_electron*0; % Neglect those ?
    delta_mu_H2_anode  = -R*temp*log(1-J_H2_anode*t_GDL/(D_H2H2O*c*moleFraction(gas_anode,'H2')));
    delta_mu_H2O_anode = R*temp*log(1-J_H2O_anode*t_GDL/(D_H2H2O*c*moleFraction(gas_anode,'H2O')));

    % Chemical potential of O2 ions (anode)
    R_eq = j/(2*F);
    mu_O2ion_anode = mu_O2ion_YSZ_anode_eq + delta_mu_H2_anode + R*temp*log(v/R_eq +...
        exp((delta_mu_H2O_anode + 2*delta_mu_e_anode)/(R*temp) ));

    % YSZ
    J_O2ion_YSZ = v;
    cond_prime_YSZ = cond_YSZ/(F^2);
    delta_mu_O2ion_YSZ = J_O2ion_YSZ*t_e/cond_prime_YSZ;
    mu_O2ion_cathode = mu_O2ion_anode + delta_mu_O2ion_YSZ;

    % Cathode side
    J_02_cathode = 0.5*v;
    x_O2_gas_cathode = moleFraction(gas_cathode,'O2');
    x_O2_cathode = 1-(1-x_O2_gas_cathode)*exp( (J_02_cathode*t_GDL) / (c*D_O2N2));
    delta_mu_O2_cathode = -R*temp*log(x_O2_cathode/x_O2_gas_cathode);
    mu_O2_cathode = mu_cathode_eq(iO2) - delta_mu_O2_cathode;

    % Cathode electron potential
    mu_electron_cathode = mu_electron_cathode_eq + 0.25*delta_mu_O2_cathode + 0.5*R*temp*...
        log(v/R_eq + exp( (mu_O2ion_cathode-mu_O2ion_YSZ_cathode_eq) / (R*temp) ));

    % Finally going to cathode side terminal
    delta_mu_e_cathode = J_electron * 0; % Neglect those ?
    mu_electron_cathode_term = mu_electron_cathode + delta_mu_e_cathode;

    delta_phi(i) = -1/F * (-mu_electron_anode_eq + mu_electron_cathode_term); %V
end

%% Plots
clf
plot(current_densities,delta_phi,'-k')
legend('Elec. Potential (V)')
xlabel("Current density (kA/m²)")
plotfixer
