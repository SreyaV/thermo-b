% Project 9: SOFC "button cell" calculations, Parts 2 and 3.
% C.F. Edwards, 2-28-09

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
e = 1.60218e-19;   % Coul/unit-charge
N_A = 6.02214e23;    % particles/mole
R   = 8.31447;       % J/mol-K
h   = 6.62607e-34;   % J-s
F   = e*N_A;         % Coul/mole-of-unit-charge and (J/mol)/(eV/particle)
k_B = R/N_A;         % J/particle-K

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
muH2a_eq  = muanode(H2);
muH2Oa_eq = muanode(H2O);
muO2c_eq  = mucathode(O2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with an equilibrium pass (zero current/reaction rate).

% Use the electrochemical potential of an electron in the anode as a
% reference value.  Set this to zero.
muEa_eq = 0;
% The oxygen ion electrochemical potential is determined by the 
% affinity being zero for the anode reaction.
muOa_eq = muH2Oa_eq + 2*muEa_eq - muH2a_eq; 
% The oxygen ion has the same electrochemical potential at the cathode when
% in equilibrium.
muOc_eq = muOa_eq;
% The electrochemical potential of an electron at the cathode is given by
% zero affinity for the cathode reaction.
muEc_eq = (1/2)*muOc_eq -(1/4)*muO2c_eq;
% The electrical potential of the cathode wrt the anode is the same as the
% electrochemical potential difference of the anode wrt to the cathode
% divided by zF.
phi_eq = (muEc_eq - muEa_eq)/(-1*F)
% Get the overall Gibbs function difference (affinity) imposed on the cell:
Affinity_eq = (muH2a_eq + (1/2)*muO2c_eq - muH2Oa_eq)
% Save 1000 C value for use in normalization later on.
Affinity_1000 = Affinity_eq

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
L_Nia  = 0;
L_Nic  = 0;
L = [L_Nia L_GDLa L_GDLa L_YSZ L_GDLc L_Nic];

% Set the gas conductivities in terms of the binary diffusivities.
concentration  = Pcell/(R*Tcell);                           % mol/m3
D_O2_N2  = Binary_Diffusivity_ijTP(iO2,iN2,Tcell,Pcell)     % m2/s
D_H2_H2O = Binary_Diffusivity_ijTP(iH2,iH2O,Tcell,Pcell)
% Save 1000 C values for use in normalization later.
D_O2_N2_1000 = D_O2_N2
D_H2_H2O_1000 = D_H2_H2O

% Work out the temperature dependence of the properties.
steps = 200;

% Ionic conductivity of YSZ via Arrhenius.
kappa_YSZo = 15;            % S/m
T_YSZo     = 1000+273.15;   % K
Ea_YSZ     = 1;             % eV/particle
Ea_YSZ     = Ea_YSZ*F;      % J/mol

% Exchange current density of the cathode.
ioco   = 1000;                      % C/s/m2
T_ioco = 1000+273.15;               % K
Ea_ioc = 100e3;                     % J/mol

i = 1;
TmaxC = 1300;
TminC = 600;
Pcell = 1e5;
for TC=TminC:10:TmaxC
    TC
    % Set the temperature.
    TCi(i) = TC;
    Tcell = TC+273.15;
    
    % Temperature-dependent conductivity for YSZ.
    kappa_YSZ = kappa_YSZo*exp(-(Ea_YSZ/R)*(1/Tcell - 1/T_YSZo));
    kappa_YSZi(i) = kappa_YSZ;
    kappa_O = kappa_YSZ/(2*F)^2;
    % Set the electron specific conductivity of the nickel.  It is too
    % small to warrant adjusting with temperature.
    kappa_Ni  = 2e6;                % S/m or (A/m2)/(V/m)
    kappa_E = kappa_Ni/(-1*F)^2;    % (mol-E/m2/s)/(J/mol-E/m)
    % Set the gas conductivities in terms of the binary diffusivities.
    concentration  = Pcell/(R*Tcell);    % mol/m3
    D_O2_N2  = Binary_Diffusivity_ijTP(iO2,iN2,Tcell,Pcell);
    D_H2_H2O = Binary_Diffusivity_ijTP(iH2,iH2O,Tcell,Pcell);
    D_O2_N2i(i)  = D_O2_N2;
    D_H2_H2Oi(i) = D_H2_H2O;  
    % Store for passing to routine.
    K = [kappa_E kappa_O concentration D_O2_N2 D_H2_H2O];
    
    % Set the gas object to the anode state:
    set(gas,'T',Tcell,'P',Pcell,'X',xanode);
    % Get the chemical potentials:
    muanode = chemPotentials(gas);      % J/kmol
    muanode = muanode/1000;             % J/mol
    % Set the gas object to the cathode state:
    set(gas,'T',Tcell,'P',Pcell,'X',xcathode);
    % Get the chemical potentials:
    mucathode = chemPotentials(gas);    % J/kmol
    mucathode = mucathode/1000;         % J/mol
    % Save the needed  chemical potentials in their own variable names.
    muEa_eq   = 0;
    muH2a_eq  = muanode(H2);
    muH2Oa_eq = muanode(H2O);
    muO_eq    = muH2Oa_eq + 2*muEa_eq - muH2a_eq;
    muO2c_eq  = mucathode(O2);
    muEc_eq   = (1/2)*muO_eq -(1/4)*muO2c_eq;
    mu_eq = [muEa_eq muH2a_eq muH2Oa_eq muO_eq muO2c_eq muEc_eq];
    % Save the needed mole fractions.
    xH2a_eq  = xanode(H2);
    xH20a_eq = xanode(H2O);
    xO2c_eq  = xcathode(O2);
    x_eq = [xH2a_eq xH20a_eq xO2c_eq];

    % Equilibrium potential and affinity.
    phi_eq = (muEa_eq - muEc_eq)/F;
    Affinityi(i) = (muH2a_eq + (1/2)*muO2c_eq - muH2Oa_eq);
    
    % Exchange current density.
    ioc = ioco*exp(-(Ea_ioc/R)*(1/Tcell - 1/T_ioco));
    ioa = Anode_Cathode_Ratio*ioc;
    ioci(i) = ioc;
    ioai(i) = ioa;
    
    % Find the peak power for this temperature.
    imax = 100*1000;    % current density, A/m2
    steps = 6000;
    vmax = imax/(2*F);  % reaction velocity density, rxn/s/m2
    vmin = 0;
    dv = (2*vmax-vmin)/steps;
    Pmax = 0;
    for v=vmin:dv:2*vmax
        v
        [phi mu] = SOFC_Element_icTKL(2*F*v,x_eq,mu_eq,Tcell,K,L,ioa,ioc);
        if(phi == 0)
            break
        end
        % Save conditions at the peak power point.
        Power = phi*2*F*v;
        Pmax = max(Pmax,Power);
        if(Pmax == Power)
            iatmax = 2*F*v;
            phiatmax = phi;
        end
    end
    Pmaxi(i)      = Pmax;
    phiatmaxi(i)  = phiatmax;
    iatmaxi(i)    = iatmax;
    etaIatmaxi(i) = Pmax/((iatmax/(2*F))*120000*2.016);

    if(TC == 1000)
        Pmax_1000 = Pmax;
        phiatPmax_1000 = phiatmax;
        iiatPmax_1000  = iatmax;
    end
    
    i = i+1;
end

figure(1)
clf
hold on
plot(TCi,Affinityi/Affinity_1000,'-','Color',[0 .5 0])
plot(TCi,D_H2_H2Oi/D_H2_H2O_1000,'--','Color',[.5 .25 0])
plot(TCi,D_O2_N2i/D_O2_N2_1000,'--','Color',[.25 .5 .0])
plot(TCi,ioci/ioco,'m--')
plot(TCi,kappa_YSZi/kappa_YSZo,'b--')
plot(TCi,Pmaxi/Pmax_1000,'k')
plot(TCi,iatmaxi/iiatPmax_1000,'g-')
plot(TCi,phiatmaxi/phiatPmax_1000,'r')
hold off
xlabel('T (\circC)')
legend('Overall Affinity',...
    'D H2/H2O',...
    'D O2/N2',...
    'i_o, cathode',...
    'O^= Conductivity',...
    'Peak Power',...
    'i at Peak Power',...
    '\Delta\phi at Peak Power')
axis([TminC TmaxC 0 3])
plotfixer


figure(2)
clf
hold on
plot(TCi,Affinityi/Affinity_1000,'-','Color',[0 .5 0])
plot(TCi,D_H2_H2Oi/D_H2_H2O_1000,'--','Color',[.5 .25 0])
plot(TCi,D_O2_N2i/D_O2_N2_1000,'--','Color',[.25 .5 .0])
plot(TCi,ioci/ioco,'m--')
plot(TCi,kappa_YSZi/kappa_YSZo,'b--')
hold off
xlabel('T (\circC)')
legend('Overall Affinity',...
    'D H2/H2O',...
    'D O2/N2',...
    'i_o, cathode',...
    'O^= Conductivity')
axis([TminC TmaxC 0 3])
plotfixer

% Do a series of V-i curves at different Ts.  Uses equilibrium info from
% above, so don't separate these.

clear phii Pi ii
j = 1;
TmaxC = 1100;
TminC = 700;
for TC=TminC:100:TmaxC
    TC
    Tcell = TC+273.15;
    
    % Temperature-dependent conductivity for YSZ.
    kappa_YSZ = kappa_YSZo*exp(-(Ea_YSZ/R)*(1/Tcell - 1/T_YSZo));
    kappa_O = kappa_YSZ/(2*F)^2;
    % Set the electron specific conductivity of the nickel.  It is too
    % small to warrant adjusting with temperature.
    kappa_Ni  = 2e6;                % S/m or (A/m2)/(V/m)
    kappa_E = kappa_Ni/(-1*F)^2;    % (mol-E/m2/s)/(J/mol-E/m)
    % Set the gas conductivities in terms of the binary diffusivities.
    concentration  = Pcell/(R*Tcell);    % mol/m3
    D_O2_N2  = Binary_Diffusivity_ijTP(iO2,iN2,Tcell,Pcell);
    D_H2_H2O = Binary_Diffusivity_ijTP(iH2,iH2O,Tcell,Pcell);
    % Store for passing to routine.
    K = [kappa_E kappa_O concentration D_O2_N2 D_H2_H2O];
    
    % Set the gas object to the anode state:
    set(gas,'T',Tcell,'P',Pcell,'X',xanode);
    % Get the chemical potentials:
    muanode = chemPotentials(gas);      % J/kmol
    muanode = muanode/1000;             % J/mol
    % Set the gas object to the cathode state:
    set(gas,'T',Tcell,'P',Pcell,'X',xcathode);
    % Get the chemical potentials:
    mucathode = chemPotentials(gas);    % J/kmol
    mucathode = mucathode/1000;         % J/mol
    % Save the needed  chemical potentials in their own variable names.
    muEa_eq   = 0;                      % Reference potential for system
    muH2a_eq  = muanode(H2);
    muH2Oa_eq = muanode(H2O);
    muO_eq    = muH2Oa_eq + 2*muEa_eq - muH2a_eq;
    muO2c_eq  = mucathode(O2);
    muEc_eq   = (1/2)*muO_eq -(1/4)*muO2c_eq;
    mu_eq = [muEa_eq muH2a_eq muH2Oa_eq muO_eq muO2c_eq muEc_eq];
    % Save the needed mole fractions.
    xH2a_eq  = xanode(H2);
    xH20a_eq = xanode(H2O);
    xO2c_eq  = xcathode(O2);
    x_eq = [xH2a_eq xH20a_eq xO2c_eq];

    % Equilibrium potential and affinity.
    phi_eq = (muEa_eq - muEc_eq)/F;
    
    % Exchange current density.
    ioc = ioco*exp(-(Ea_ioc/R)*(1/Tcell - 1/T_ioco));
    ioa = Anode_Cathode_Ratio*ioc;
    
    % Sweep out V-i and P-i curves
    imax = 100*1000;
    vmax = imax/2/F;
    vmin = 0;
    steps = 1000;
    dv = (vmax-vmin)/steps;
    i = 1;
    for v=vmin:dv:vmax 
        [phi mu xac] = SOFC_Element_icTKL(2*F*v,x_eq,mu_eq,Tcell,K,L,ioa,ioc);
        if(phi == 0)
            ii(j,i)   = ii(j,i-1);
            phii(j,i) = phii(j,i-1);
            Pi(j,i)   = Pi(j,i-1);
            etai(j,i) = etai(j,i-1);
            xH2ai(j,i)  = xH2ai(j,i-1);
            xH2Oai(j,i) = xH2Oai(j,i-1);
            xO2ci(j,i)  = xO2ci(j,i-1);
        else
            ii(j,i)   = 2*F*v;
            phii(j,i) = phi;
            Pi(j,i)   = phi*2*F*v;
            etai(j,i) = phi*2*F/(120000*2.016);
            xH2ai(j,i)  = xac(1);
            xH2Oai(j,i) = xac(2);
            xO2ci(j,i)  = xac(3);
        end
        i = i+1;
    end
    TCj(j) = TC;
    [Pmax(j) iatmax(j)] = max(Pi(j,:));
    j = j+1;
end
jmax = j-1;

figure(3)
clf
hold on
for j=1:1:jmax
    plot(ii(j,:)/1000,Pi(j,:)/40000,'b')
    plot(ii(j,:)/1000,phii(j,:),'k')
%     plot(ii(j,:)/1000,etai(j,:),'Color',[0 .75 0])
    text(ii(j,iatmax(j))/1000,Pmax(j)/40000,...
        [sprintf('%d',TCj(j)) '\circC'],'Color','b')
end
hold off
xlabel('Current Density (kA/m^2)')
scale = axis;
% axis([0 scale(2) 0 1.2])
legend(...
    'Power Density (40kW/m2)',...
    'Output Potential (V)')
plotfixer

figure(4)
clf
hold on
for j=1:1:jmax
    plot(ii(j,:)/1000,Pi(j,:)/40000,'b')
%     plot(ii(j,:)/1000,phii(j,:),'k')
    plot(ii(j,:)/1000,etai(j,:),'Color',[0 .5 0])
    text(ii(j,iatmax(j))/1000,Pmax(j)/40000,...
        [sprintf('%d',TCj(j)) '\circC'],'Color','b')
end
hold off
xlabel('Current Density (kA/m^2)')
scale = axis;
axis([0 scale(2) 0 1.1])
legend(...
    'Power Density (40kW/m2)',...
    'Efficiency (LHV)')
plotfixer

figure(5)
clf
hold on
for j=1:1:jmax
    plot(ii(j,:)/1000,xH2ai(j,:),'r')
    plot(ii(j,:)/1000,xO2ci(j,:),'Color',[0 .75 0])
    plot(ii(j,:)/1000,xH2Oai(j,:),'b')
end
hold off
xlabel('Current Density (kA/m^2)')
ylabel('Mole Fraction')
scale = axis;
% axis([0 scale(2) 0 1.2])
legend(...
    'Anode Hydrogen',...
    'Cathode Oxygen',...
    'Anode Water')
plotfixer
