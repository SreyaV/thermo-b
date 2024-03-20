close all;clear;clc;
format long e;
%%

% Global Variables
global e N_A R h F k_B
e = 1.60218e-19;    % Coul/unit-charge
N_A = 6.02214e23;   % particles/mole
R   = 8.31447;      % J/mol-K
h   = 6.62607e-34;  % J-s
F   = e*N_A;        % Coul/mole-of-unit-charge and (J/mol)/(eV/particle)
k_B = R/N_A;        % J/particle-K


%% make plots

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


%% 
% Evaluate first button cell differential element
% walk it to the current
phi_eq = SOFC_Element_icTKL(0,x_eq,mu_eq,Tcell,K,L,ioa,ioc);


%%
nVolts = 25;
voltageVec = linspace(phi_eq,0.001,nVolts);
currSteps = 50;
for viter = 1:1:length(voltageVec)
    [currentVals(viter),power(viter),heat_transfer(viter),fuel_utilization(viter),a,LHV_H2_unused(viter),LHV_eff_input(viter),LHV_eff_consumed(viter)] = ...
        channel_LHV(voltageVec(viter),currSteps);

    while isnan(currentVals(viter))
        currSteps = currSteps + 10;
        currentVals(viter) = channel(voltageVec(viter),currSteps);
    end
end

% second thingy
[~,maxiter] = max(power);
maxP_voltage = voltageVec(maxiter);
currSteps = 100;
[maxP_current,maxP,maxP_heat_tranasfer,maxP_fuel_utilization,maxP_alongChannel] = ...
        channel(maxP_voltage,currSteps);


%% make plots
% figure(1)
% plot(currentVals,voltageVec,'k-','LineWidth',2);
% hold on;
% plot(currentVals,power/100,'b-','LineWidth',2);
% plot(currentVals,heat_transfer/100,'r-','LineWidth',2);
% plot(currentVals,fuel_utilization,'-','Color',[0 0.5 0],'LineWidth',2);
% 
% ylim([0 1.75])
% 
% xlabel('Cell current (A)');
% 
% lgd = legend('Electric potential (V)','Power output/100 (W)','Heat transfer/100 (W)','Fuel utilization');
% 
% set(lgd,'Location','NorthWest','Box','Off');
% 
% set(gca,'FontSize',26,'FontWeight','Bold');
% set(gcf,'Color','white','Position',[680 172 975 706]);
% 
% figure(2);clf;
% plot(maxP_alongChannel.distance_along_channel,...
%     maxP_alongChannel.current_density_array/1000,'k-','LineWidth',2);
% hold on;
% plot(maxP_alongChannel.distance_along_channel,...
%     maxP_alongChannel.electrical_power_density/1000,'b-','LineWidth',2);
% plot(maxP_alongChannel.distance_along_channel,...
%     abs(maxP_alongChannel.heat_flux_array/1000),'r-','LineWidth',2);
% 
% 
% set(gca,'FontSize',26,'FontWeight','Bold');
% set(gcf,'Color','white','Position',[680 172 975 706]);
% 
% 
% 
% figure(3);clf;
% hold on;
% plot(maxP_alongChannel.distance_along_channel, maxP_alongChannel.hydrogen_mole_fractions,'k-','LineWidth',2)
% plot(maxP_alongChannel.distance_along_channel, maxP_alongChannel.oxygen_mole_fractions,'b-','LineWidth',2)
% plot(maxP_alongChannel.distance_along_channel, maxP_alongChannel.water_mole_fractions,'r-','LineWidth',2)
% plot(maxP_alongChannel.distance_along_channel, maxP_alongChannel.equ_electric_potential,'m-','LineWidth',2)
% plot(maxP_alongChannel.distance_along_channel, maxP_alongChannel.actual_electric_potential, 'k--','Linewidth',2)
% 
% xlabel("Distance Along Channel (m)")
% title(['1000 C, 3 bar, ' num2str(round(maxP_voltage,3)) ' V'])
% 
% legend("xH2 (anode)", "xO2 (cathode)", "xH2O (anode)", "\Delta \phi equil. (V)", "\Delta \phi cell (V)")

%%
figure(2)
% LHV_H2_unused = molar_flow_rate_H2 * M_H2 * LHV_hyd /200;
% LHV_eff_input = power / (initial_H2*M_H2*LHV_hyd);
% LHV_eff_consumed = power/ ((-molar_flow_rate_H2 + initial_H2) *LHV_hyd*M_H2);
plot(currentVals,voltageVec,'-k')
hold on
plot(currentVals,LHV_H2_unused,'-r')
plot(currentVals,LHV_eff_input,'-b')
plot(currentVals,LHV_eff_consumed,'-g')
legend("Elec. Potential(V)","H2 HV/200 (W)","LHV efficiency (Hyd. in)","LHV efficiency (Hyd. used)")
xlabel("Cell Current (A)")
hold off
plotfixer










