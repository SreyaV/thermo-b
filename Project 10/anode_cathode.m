function [x mu] = anode_cathode(x_H2, x_H2O, x_O2, x_N2)
% Takes mole fractions of the anode and cathode gases, returns
% composition and chemical potential arrays

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

end
