function [phi mu xac] = SOFC_Element_icTKL(i,x_eq,mu_eq,Tcell,K,L,ioa,ioc)
% Return the electrical potential difference between the cathode and anode
% (V) and the array of electrochemical potentials (J/mol) for a
% differential area element of a hydrogen-water-air SOFC operating at
% current density i (A/m2) with channel (and therefore equilibrium) 
% chemical potentials in array mu_eq (J/mol), conductivities K 
% (mol-i/m2/s)/(J/mol-i), path lengths L (m), and anode and cathode
% exchange current densities ioa and ioc (A/m2/s).

global R F

% Set a logic variable to control output comments.
verbose = false;

% Return zeros for the mole fractions if not otherwise found.
xac = [0 0 0];

% Unpack the array of equilibrium (channel) mole fractions.
xH2a_eq  = x_eq(1);
xH2Oa_eq = x_eq(2);
xO2c_eq  = x_eq(3);

% Unpack the array of equilibrium (channel) e-c potentials.
muEa_eq   = mu_eq(1);
muH2a_eq  = mu_eq(2);
muH2Oa_eq = mu_eq(3);
muO_eq    = mu_eq(4);
muO2c_eq  = mu_eq(5);
muEc_eq   = mu_eq(6);
muOa_eq   = muO_eq;
muOc_eq   = muO_eq;

% Unpack the conductivities and transport path lengths.
kappa_E  = K(1);
kappa_O  = K(2);
Conc     = K(3);
D_O2_N2  = K(4);
D_H2_H2O = K(5);

L_Ea  = L(1);
L_H2  = L(2);
L_H2O = L(3);
L_O   = L(4);
L_O2  = L(5);
L_Ec  = L(6);

% Set the equilibrium reaction rate densities (based on the current
% exchange densities) and the net reaction rate (velocity) based on i.
Roa = ioa/(2*F);    % mol-rxn/s/m2
Roc = ioc/(2*F);    % mol-rxn/s/m2
v   = i/(2*F);      % mol-rxn/s/m2

% Must draw H2 from gas supply.
JH2a = v;   % Flux, mol-i/m2-s
% Change in mole fraction across GDL.
Delta_xH2a  = JH2a*L_H2/(Conc*D_H2_H2O);
xH2a = xH2a_eq - Delta_xH2a;
if(xH2a <= 0)
    if(verbose)
        disp('All H2 consumed...')
    end
    phi = 0;
    mu = zeros(1,9);
    return
end
% Change in chemical potential across GDL.
Delta_muH2a = -R*Tcell*log(1 - JH2a*L_H2/(xH2a_eq*Conc*D_H2_H2O));
muH2a  = muH2a_eq - Delta_muH2a;

% Must drive H2O to gas supply.
JH2Oa = v;  % Flux, mol-i/m2-s
% Change in mole fraction across GDL.
Delta_xH2Oa  = JH2Oa*L_H2O/(Conc*D_H2_H2O);
xH2Oa = xH2Oa_eq + Delta_xH2Oa;
if(xH2Oa >= 1)
    if(verbose)
        disp('Anode H2O hit unity...')
    end
    phi = 0;
    mu = zeros(1,9);
    return
end
% Change in chemical potential across GDL.
Delta_muH2Oa = R*Tcell*log(1 + JH2Oa*L_H2O/(xH2Oa_eq*Conc*D_H2_H2O) );
muH2Oa = muH2Oa_eq + Delta_muH2Oa;

% Choose the anode electron at the terminal as reference.
muEat  = muEa_eq;   % Use equilibrium reference value.
% Must drive electrons to the terminal.
JEa = 2*v;  % Flux, mol-i/m2-s
Delta_muEa = JEa*L_Ea/kappa_E;
muEa   = muEat + Delta_muEa;

% The oxygen ion e-c potential comes from anode reaction.
muOa   = muOa_eq + Delta_muH2a + (R*Tcell)*log(...
    v/Roa + exp((Delta_muH2Oa + 2*Delta_muEa)/(R*Tcell))...
    );

% Must drive oxygen ion across the YSZ.
JOac = v;   % Flux, mol-i/m2-s
Delta_muOac = JOac*L_O/kappa_O;
muOc   = muOa + Delta_muOac;

% Must draw O2 from gas supply.
JO2c = (1/2)*v; % Flux, mol-i/m2-s
% Mole fraction at cathode.
xO2c = 1 - (1-xO2c_eq)*exp(JO2c*L_O2/(Conc*D_O2_N2));
if(xO2c <= 0)
    if(verbose)
        disp('All O2 consumed...')
    end
    phi = 0;
    mu = zeros(1,9);
    return
end
% Change in chemical potential across GDL.
Delta_muO2c = -R*Tcell*log(xO2c/xO2c_eq);
muO2c  = muO2c_eq - Delta_muO2c;

% The electron e-c potential comes from cathode reaction.
muEc   = muEc_eq + (1/4)*Delta_muO2c + (R*Tcell/2)*log(...
    v/Roc + exp((muOc - muOc_eq)/(R*Tcell))...
    );

% Must draw the electron from the terminal.
JEc = 2*v;  % Flux, mol-i/m2-s
Delta_muEc = JEc*L_Ec/kappa_E;
muEct  = muEc + Delta_muEc;

% Return the electrical potential across the terminals and a mu array.
phi = -(muEct-muEat)/F;
mu = [muEat muEa muH2a muH2Oa muOa muOc muO2c muEc muEct];

% Can return the mole fractions at the electrodes too if you want.
% This is added to answer Simon's questions about negative potentials.
% Might need to remove this for compatibility reasons next year.  Not sure...
xac = [xH2a xH2Oa xO2c];
