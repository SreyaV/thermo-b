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

delta_phi_eq = 0.5/F * ( mu_anode_eq(iH2) + 0.5 * mu_cathode_eq(iO2) -mu_anode_eq(iH2O))


%% Pass 2

disp("Pass 2: Out of Equilibrium ...")
dj = 0.01;
jmax = 55;
current_densities = 1e-9:dj:jmax;

for i=1:length(current_densities)
    j = current_densities(i)*1000; % Current density A/m²
    v = (1/2/F)*j; % Reaction velocity per area : mol/s/m²

    % Fluxes per area: anode
    J_H2_anode  = v;
    J_H2O_anode = v;
    J_electron  = 2*v;

    % Anode activation losses
    c = P/ (R*temp); % Total concentration of species
    delta_mu_e_anode = J_electron*0; % Neglect those ?
    
    % Anode GDL losses
    x_H2_gas_anode = moleFraction(gas_anode,'H2');
    H2_anode(i) = x_H2_gas_anode*( 1 - (J_H2_anode*t_GDL)/(x_H2_gas_anode*c*D_H2H2O));
    x_H2O_gas_anode = moleFraction(gas_anode,'H2O');
    H2O_anode(i) = x_H2O_gas_anode*( 1 + (J_H2O_anode*t_GDL)/(x_H2O_gas_anode*c*D_H2H2O));

    delta_mu_H2_anode  = -R*temp*log(H2_anode(i) / x_H2_gas_anode);
    delta_mu_H2O_anode = R*temp*log(H2O_anode(i) / x_H2O_gas_anode);
    mu_H2_anode = mu_anode_eq(iH2) - delta_mu_H2_anode;
    mu_H2O_anode = mu_anode_eq(iH2O) + delta_mu_H2O_anode;
    GDL_loss_anode(i) = (delta_mu_H2_anode + delta_mu_H2O_anode)/F;

    % Chemical potential of O2 ions (anode)
    R_anode_eq = J_anode/(2*F);
    mu_O2ion_anode = mu_O2ion_YSZ_anode_eq + delta_mu_H2_anode + R*temp*log(v/R_anode_eq +...
        exp((delta_mu_H2O_anode + 2*delta_mu_e_anode)/(R*temp) ));
    anode_loss(i) = (mu_O2ion_anode - (mu_H2O_anode - mu_H2_anode))/F;

    % YSZ loss
    J_O2ion_YSZ = v;
    cond_prime_YSZ = cond_YSZ/((2*F)^2);
    delta_mu_O2ion_YSZ = J_O2ion_YSZ*t_e/cond_prime_YSZ;
    mu_O2ion_cathode = mu_O2ion_anode + delta_mu_O2ion_YSZ;
    ohmic_loss(i) = delta_mu_O2ion_YSZ/F;

    % Cathode GDL losses
    J_02_cathode = 0.5*v;
    x_O2_gas_cathode = moleFraction(gas_cathode,'O2');
    x_O2_cathode = 1-(1-x_O2_gas_cathode)*exp( (J_02_cathode*t_GDL) / (c*D_O2N2));

    O2_cathode(i) = x_O2_cathode;
    delta_mu_O2_cathode = -R*temp*log(x_O2_cathode/x_O2_gas_cathode);
    mu_O2_cathode = mu_cathode_eq(iO2) - delta_mu_O2_cathode;
    GDL_loss_cathode(i) = 0.5*delta_mu_O2_cathode/F;

    % Cathode electron potential
    R_cathode_eq = J_cathode/(2*F);
    mu_electron_cathode = mu_electron_cathode_eq + 0.25*delta_mu_O2_cathode + 0.5*R*temp*...
        log(v/R_cathode_eq + exp( (mu_O2ion_cathode-mu_O2ion_YSZ_cathode_eq) / (R*temp) ));
    cathode_loss(i) = 2*(mu_electron_cathode - (0.5*mu_O2ion_cathode-0.25*mu_O2_cathode) )/F;

    % Cathode Activation losses
    delta_mu_e_cathode = J_electron * 0; % Neglect those ?
    mu_electron_cathode_term = mu_electron_cathode + delta_mu_e_cathode;
    
    % Elec Potential (V)
    delta_phi(i) = -1/F * (-mu_electron_anode_eq + mu_electron_cathode_term); %V

    % Power density
    power(i) = delta_phi(i)*j; % (W/m²)

    GDL_loss(i) = GDL_loss_cathode(i) + GDL_loss_anode(i);
    
    % Check that we didn't reach limit
    if imag(delta_phi(i))>1e-6
        break;
    end
end

% Re-size current_densities
current_densities = current_densities(1:length(delta_phi));

%% Plots: Potential & co
% Take all arrays up to end-1 to remove the 1st complex number

figure(1)
plot(current_densities,delta_phi,'-k')
hold on
plot(current_densities(1:end-1),power(1:end-1)/max(power),'-',color='cyan')
plot(current_densities(1:end-1),H2_anode(1:end-1),'-r')
plot(current_densities(1:end-1),H2O_anode(1:end-1),'-b')
plot(current_densities(1:end-1),O2_cathode(1:end-1),'-g')
axis([0 60 0 1.2])
text(2,1.1,"Max Power: "+ string(round(max(power)/1e3,2))+" kW/m²")

legend('Elec. Potential (V)','Power (normalized)', 'Anode Hydrogen', 'Anode Water', 'Cathode Oxygen')
xlabel("Current density (kA/m²)")
hold off
plotfixer

%% Losses plot

figure(2)
plot(current_densities(1:end-1),GDL_loss(1:end-1),'-g')
hold on
plot(current_densities(1:end-1),ohmic_loss(1:end-1),'-b')
plot(current_densities(1:end-1),anode_loss(1:end-1),'-r')
plot(current_densities(1:end-1),cathode_loss(1:end-1),'-k')

legend('GDL loss','Ohmic loss (YSZ)','Anode Loss','Cathode Loss')
xlabel("Current density (kA/m²)")
ylabel('Electrochem. Potential (eV/rxn)')
hold off
plotfixer

%% Affinity plot

for i=1:length(current_densities)
    affinity0(i) = delta_phi_eq*2; % factor 2 because 2 e-
end
affinity1 = affinity0 - anode_loss;
affinity2 = affinity1 - cathode_loss;
affinity3 = affinity2 - ohmic_loss;
affinity4 = affinity3 - GDL_loss;

figure(3)
plot(current_densities(1:end-1),affinity0(1:end-1),'-b')
hold on
plot(current_densities(1:end-1),affinity1(1:end-1),'-r')
plot(current_densities(1:end-1),affinity2(1:end-1),'-r')
plot(current_densities(1:end-1),affinity3(1:end-1),'-g')
plot(current_densities(1:end-1),affinity4(1:end-1),'-k')
plot(current_densities(1:end-1),delta_phi(1:end-1)*2,'-k') %Output potential x2 (to check)

legend('Affinity (reactants)','Anode Loss','Cathode Loss','Ohmic Loss','GDL Loss')
xlabel("Current density (kA/m²)")
ylabel('Electrochem. Potential (eV/rxn)')
hold off
plotfixer

