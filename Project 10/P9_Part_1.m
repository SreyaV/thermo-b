% Project 9: SOFC "button cell" calculations, Part 1.
% C.F. Edwards, 2-28-10

clear all
format compact
fprintf('\n************************************************************\n')

% Species designation in function Binary_Diffusivity_ijTP:
iN2 = 1; iO2 = 2; iAr = 3; iCO2 = 4; iH2O = 5; iCO = 6; iH2 = 7; iCH4 = 8;

% Some constants that are hard to remember:
% 1 J = 1 Coul-Volt
% 1 eV = 1.602e-19 J
% F = e*N_A = 96485.3 Coul/mole-of-unit-charge

global e N_A R h F k_B
e = 1.60218e-19;    % Coul/unit-charge
N_A = 6.02214e23;   % particles/mole
R   = 8.31447;      % J/mol-K
h   = 6.62607e-34;  % J-s
F   = e*N_A;        % Coul/mole-of-unit-charge and (J/mol)/(eV/particle)
k_B = R/N_A;        % J/particle-K

% Get a Cantera ideal gas object.  
% Type this for a list of properties:  methods solution -full
gas = Solution('GRI30.yaml');
H2  = speciesIndex(gas,'H2');
H2O = speciesIndex(gas,'H2O');
O2  = speciesIndex(gas,'O2');
N2  = speciesIndex(gas,'N2');
Nsp = nSpecies(gas);

% Set the temperature and pressure of the (isothermal, iobaric) cell:
Tcell = 1000+273.15
Pcell = 1e5

% Set the composition and pressure at the anode:
xanode = zeros(Nsp,1);
xanode(H2O) = 0.03;
xanode(H2)  = 1-xanode(H2O);
% Set the gas object to the anode state:
set(gas,'Temperature',Tcell,'Pressure',Pcell,'MoleFractions',xanode);
% Get the chemical potentials:
muanode = chemPotentials(gas);      % J/kmol
muanode = muanode/1000;             % J/mol

% Set the composition and pressure at the cathode:
xcathode = zeros(Nsp,1);
xcathode(O2) = 0.21;
xcathode(N2) = 1-xcathode(O2);
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
muH2a_eq  = muanode(H2)
muH2Oa_eq = muanode(H2O)
muO2c_eq  = mucathode(O2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with an equilibrium pass (zero current/reaction rate).

% Use the electrochemical potential of an electron in the anode as a
% reference value.  Set this to zero.
muEa_eq = 0
% The oxygen ion electrochemical potential is determined by the 
% affinity being zero for the anode reaction: H2 + O= -> H2O + 2e-
muOa_eq = muH2Oa_eq + 2*muEa_eq - muH2a_eq
% The oxygen ion has the same electrochemical potential at the cathode when
% in equilibrium.
muOc_eq = muOa_eq
% The electrochemical potential of an electron at the cathode is given by
% zero affinity for the cathode reaction: O2 + 4e- -> 2O=
muEc_eq = (2/4)*muOc_eq -(1/4)*muO2c_eq
% The electrical potential of the cathode wrt the anode is the same as the
% electrochemical potential difference of the anode wrt to the cathode
% divided by zF.
phi_eq = (muEc_eq - muEa_eq)/(-1*F)
% Get the overall Gibbs function difference (affinity) imposed on the cell:
Affinity_eq = (muH2a_eq + (1/2)*muO2c_eq - muH2Oa_eq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now do it with a specified, finite-rate loss.  

% Set the cathode exchange current density (and thereby eq. reaction rate).
ioc = 1000;                     % (C/s)/m2
Roc = ioc/(2*F);                % (mol-rxn/s)/m2
% Set the anode kinetics wrt the cathode.
Anode_Cathode_Ratio = 100
Roa = Anode_Cathode_Ratio*Roc;
% Set the oxygen ion specific conductivity of the YSZ.
kappa_YSZ = 15;                 % S/m or (A/m2)/(V/m)
% Since two charges move with each O= ion, the species conductivity
% (product of concentration and molar diffusivity) expressed
% in terms of chemical potential as the driving force is this
% over (zF)^2.
kappa_O = kappa_YSZ/(-2*F)^2    % (mol-O/m2/s)/(J/mol-O/m)
% Set the electron specific conductivity of the nickel.
kappa_Ni = 2e6;                 % S/m or (A/m2)/(V/m)
kappa_E = kappa_Ni/(-1*F)^2     % (mol-E/m2/s)/(J/mol-E/m)

% Set the path lengths.
L_YSZ  = 50e-6;                 % m
L_GDLa = 5e-3;
L_GDLc = 5e-3;
L_Nia  = 0;                     % In case you want the metal too.
L_Nic  = 0;
% Put these in an array for convenient passing.
L = [L_Nia L_GDLa L_GDLa L_YSZ L_GDLc L_Nic];

% Set the gas conductivities in terms of the binary diffusivities.
concentration  = Pcell/(R*Tcell);                           % mol/m3
D_O2_N2  = Binary_Diffusivity_ijTP(iO2,iN2,Tcell,Pcell)     % m2/s
D_H2_H2O = Binary_Diffusivity_ijTP(iH2,iH2O,Tcell,Pcell)
% Save 1000C value for use in normalization later.
D_O2_N2_1000 = D_O2_N2
D_H2_H2O_1000 = D_H2_H2O

% Run a loop through the net reaction rate range.

% Set the range based on the exchange current density.
imax = 100*1000;    % current density, A/m2
vmax = imax/(2*F);  % reaction velocity density, rxn/s/m2
vmin = 0;
steps = 10000;
dv = (vmax-vmin)/steps;
i = 1;
for v=vmin:dv:vmax
    % Must draw H2 from gas supply.  Rate required is one H2 per rxn.
    JH2a = v;   % Flux, mol-i/m2-s
    % Change in mole fraction across GDL for this flux rate.
    Delta_xH2a  = JH2a*L_GDLa/(concentration*D_H2_H2O);
    xH2a = xH2a_eq - Delta_xH2a;
    % Look to see if you are out of H2.
    if(xH2a <= 0)
        disp('All H2 consumed...')
        break
    end
    % Change in chemical potential across GDL.
    Delta_muH2a = -R*Tcell*log(1 - JH2a*L_GDLa/(xH2a_eq*concentration*D_H2_H2O));
    muH2a  = muH2a_eq - Delta_muH2a;
    
    % Must drive H2O to gas supply.  Rate is one H2O per rxn.
    JH2Oa = v;  % Flux, mol-i/m2-s
    % Change in mole fraction across GDL.
    Delta_xH2Oa  = JH2Oa*L_GDLa/(concentration*D_H2_H2O);
    xH2Oa = xH2Oa_eq + Delta_xH2Oa;
    % Look to see if you are saturated with H2O.
    if(xH2Oa >= 1)
        disp('Anode H2O hit unity...')
        break
    end
    % Change in chemical potential across GDL.
    Delta_muH2Oa = R*Tcell*log(1 + JH2Oa*L_GDLa/(xH2Oa_eq*concentration*D_H2_H2O) );
    muH2Oa = muH2Oa_eq + Delta_muH2Oa;
    
    % Choose the anode electron at the terminal as reference.
    muEat  = 0;
    % Must drive electrons to the terminal.  Need two per rxn.
    JEa = 2*v;  % Flux, mol-i/m2-s
    Delta_muEa = JEa*L_Nia/kappa_E;
    muEa   = muEat + Delta_muEa;
    
    % The oxygen ion e-c potential comes from the anode reaction.
    muOa   = muOa_eq + Delta_muH2a + (R*Tcell)*log(...
        v/Roa + exp((Delta_muH2Oa + 2*Delta_muEa)/(R*Tcell))...
        );

    % Must drive oxygen ion across the YSZ.  Need one per rxn.
    JOac = v;   % Flux, mol-i/m2-s
    Delta_muOac = JOac*L_YSZ/kappa_O;
    muOc   = muOa + Delta_muOac;
    
    % Must draw O2 from gas supply.  Only need half a mole per mole rxn.
    JO2c = (1/2)*v; % Flux, mol-i/m2-s
    % Mole fraction at cathode.
    xO2c = 1 - (1-xO2c_eq)*exp(JO2c*L_GDLc/(concentration*D_O2_N2));
    % Look to see if there is O2 left.
    if(xO2c <= 0)
        disp('All O2 consumed...')
        break
    end
    % Change in chemical potential across GDL.
    Delta_muO2c = -R*Tcell*log(xO2c/xO2c_eq);
    muO2c  = muO2c_eq - Delta_muO2c;
    
    % The electron e-c potential comes from the cathode reaction.
    muEc   = muEc_eq + (1/4)*Delta_muO2c + (R*Tcell/2)*log(...
        v/Roc + exp((muOc - muOc_eq)/(R*Tcell))...
        );

    % Must draw the electron from the terminal.  Still need two per rxn.
    JEc = 2*v;  % Flux, mol-i/m2-s
    Delta_muEc = JEc*L_Nic/kappa_E;
    muEct  = muEc + Delta_muEc;
    
    % Gather together the results at this reaction velocity 
    % (current density).
    ii(i) = 2*F*v;                      % current density, A/m2
    xH2ai(i)  = xH2a;                   % Anode actual mole fraction
    xH2Oai(i) = xH2Oa;                  % Anode actual mole fraction
    xO2ci(i)  = xO2c;                   % Cathode actual mole fraction
    phiti(i) = (muEct-muEat)/(-1*F);    % cell electrical potential at terminals
    phii(i)  = (muEc-muEa)/(-1*F);      % cell electrical potential wout term loss
    poweri(i) = ii(i)*phiti(i);         % cell power density, W/m2
    % Overall affinity, J/mol-rxn
    DeltaGi(i) = muH2a + (1/2)*muO2c - muH2Oa;
    % Loss of e-c potential due to gas transport, J/mol-rxn
    GDL_loss(i)     = Delta_muH2a + (1/2)*Delta_muO2c + Delta_muH2Oa;
    % Loss of e-c potential due to ion & electron transport, J/mol-rxn
    ohmic_loss(i)   = v*L_YSZ/kappa_O + 2*v*(L_Nia+L_Nic)/kappa_E;
    % Loss of e-c potential due to activation at anode, J/mol-rxn
    anode_loss(i)   = muH2a + muOa - muH2Oa - 2*muEa;
    % Loss of e-c potential due to activation at cathode, J/mol-rxn
    cathode_loss(i) = (1/2)*muO2c + 2*muEc - muOc;
    % Actual vs. eq. reactant e-c potential at anode, J/mol-rxn
    DGraeqi(i) = muH2a + muOa - muH2a_eq - muOa_eq;
    % Actual vs. eq. product e-c potential at anode, J/mol-rxn
    DGpaeqi(i) = muH2Oa + 2*muEa - muH2Oa_eq - 2*muEa_eq;
    % Actual vs. eq. reactant e-c potential at cathode, J/mol-rxn
    DGrceqi(i) = (1/2)*muO2c + 2*muEc - (1/2)*muO2c_eq - 2*muEc_eq;
    % Actual vs. eq. product e-c potential at cathode, J/mol-rxn
    DGpceqi(i) = muOc - muOc_eq;
    % Save indivividual e-c potentials, J/mol
    muH2ai(i)  = muH2a;
    muOai(i)   = muOa;
    muH2Oai(i) = muH2Oa;
    muEai(i)   = muEa;
    muEati(i)  = muEat;
    muO2ci(i)  = muO2c;
    muOci(i)   = muOc;
    muEci(i)   = muEc;
    muEcti(i)  = muEct;
    
    i = i+1;
end

% Find the max power value and location.
[Pmax index_at_max] = max(poweri)

figure(1)
clf
hold on
plot(ii/1000,phiti,'k-')
plot(ii/1000,poweri/Pmax,'-','Color',[.5 .5 .25])
scale = axis;
axis([0 scale(2) 0 1.2])
plot(ii/1000,xH2ai,'r')
plot(ii/1000,xO2ci,'-','Color',[0 .5 0])
plot(ii/1000,xH2Oai,'b')
plot([0 scale(2)],[xH2a_eq xH2a_eq],'r--')
plot([0 scale(2)],[xO2c_eq xO2c_eq],'--','Color',[0 .5 0])
plot([0 scale(2)],[xH2Oa_eq xH2Oa_eq],'b--')
hold off
xlabel('Current Density (kA/m^2)')
legend('Elec.Potential (V)',...
       'Power/Max.Power',...
       'Anode Hydrogen',...
       'Cathode Oxygen',...
       'Anode Water')
text(2,1.1,sprintf('Max.Power: %.1f kW/m^2',Pmax/1000),'Color',[.5 .5 .25])
title(sprintf('YSZ: %.0f C, %.0f bar',Tcell-273.15,Pcell/1e5))
plotfixer

figure(2)
clf
hold on
plot([0 scale(2)],[Affinity_eq/F Affinity_eq/F],'m-')
plot(ii/1000,(Affinity_eq-GDL_loss)/F,'-','Color',[0 .5 0])
plot(ii/1000,(Affinity_eq-GDL_loss-ohmic_loss)/F,'b-')
plot(ii/1000,(Affinity_eq-GDL_loss-ohmic_loss-anode_loss)/F,'r-')
plot(ii/1000,(Affinity_eq-GDL_loss-ohmic_loss-anode_loss-cathode_loss)/F,'k-')
hold off
xlabel('Current Density (kA/m^2)')
ylabel('Electrochemical Potential (eV/rxn)')
scale = axis;
axis([0 scale(2) -.1 2.5])
legend('Reactant Affinity',...
    '  - GDL Loss',...
    '  - Ohmic Loss',...
    '  - Anode Loss',...
    '  - Cathode Loss')
title(sprintf('YSZ: %.0f C, %.0f bar',Tcell-273.15,Pcell/1e5))
plotfixer

figure(3)
clf
hold on
plot(ii/1000,GDL_loss/F,'-','Color',[0 .5 0])
plot(ii/1000,ohmic_loss/F,'b-')
plot(ii/1000,anode_loss/F,'r-')
plot(ii/1000,cathode_loss/F,'k-')
hold off
xlabel('Current Density (kA/m^2)')
ylabel('Electrochemical Potential (eV/rxn)')
scale = axis;
% axis([0 scale(2) -.1 0.5])
legend('GDL Loss','Ohmic Loss','Anode Loss','Cathode Loss')
title(sprintf('YSZ: %.0f C, %.0f bar',Tcell-273.15,Pcell/1e5))
plotfixer

figure(4)
clf
hold on
plot(ii/1000,xH2ai,'r')
plot(ii/1000,xO2ci,'-','Color',[0 .5 0])
plot(ii/1000,xH2Oai,'b')
plot([0 scale(2)],[xH2a_eq xH2a_eq],'r--')
plot([0 scale(2)],[xO2c_eq xO2c_eq],'--','Color',[0 .5 0])
plot([0 scale(2)],[xH2Oa_eq xH2Oa_eq],'b--')
hold off
xlabel('Current Density (kA/m^2)')
ylabel('Mole Fraction at Electrode')
scale = axis;
axis([0 scale(2) scale(3) scale(4)])
legend('Hydrogen','Oxygen','Water')
title(sprintf('YSZ: %.0f C, %.0f bar',Tcell-273.15,Pcell/1e5))
plotfixer

figure(5)
clf
hold on
plot(ii/1000,DGraeqi/F,'b-')
plot(ii/1000,(muH2ai-muH2a_eq)/F,'b--')
plot(ii/1000,(muOai-muOa_eq)/F,'b:')
plot(ii/1000,DGpaeqi/F,'r-')
plot(ii/1000,(muH2Oai-muH2Oa_eq)/F,'r--')
plot(ii/1000,2*(muEai-muEa_eq)/F,'r:')
hold off
xlabel('Current Density (kA/m^2)')
ylabel('Electrochemical Potential (eV/rxn)')
scale = axis;
axis([0 scale(2) scale(3) scale(4)])
legend('Anode R-Req','H2-H2eq','O-Oeq','Anode P-Peq','H2O-H2Oeq','2(E-Eeq)')
title(sprintf('YSZ: %.0f C, %.0f bar',Tcell-273.15,Pcell/1e5))
plotfixer

figure(6)
clf
hold on
plot(ii/1000,DGrceqi/F,'b-')
plot(ii/1000,(1/2)*(muO2ci-muO2c_eq)/F,'b--')
plot(ii/1000,2*(muEci-muEc_eq)/F,'b:')
plot(ii/1000,DGpceqi/F,'r-')
plot(ii/1000,(muOci-muOc_eq)/F,'r--')
hold off
xlabel('Current Density (kA/m^2)')
ylabel('Electrochemical Potential (eV/rxn)')
scale = axis;
axis([0 scale(2) scale(3) scale(4)])
legend('Cathode R-Req','(O2-O2eq)/2','2(E-Eeq)','Cathode P-Peq','O-Oeq')
title(sprintf('YSZ: %.0f C, %.0f bar',Tcell-273.15,Pcell/1e5))
plotfixer
