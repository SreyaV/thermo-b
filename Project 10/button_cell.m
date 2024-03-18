function [current_density, terminal_voltage, Pmax] = button_cell(specified_voltage, T_set, P_set)
    % Initialize constants and parameters
    clear all;  % Note: Using clear all inside a function is not standard; globals and persistent variables remain
    format compact;

    global e N_A R h F k_B;
    e = 1.60218e-19;    % Coul/unit-charge
    N_A = 6.02214e23;   % particles/mole
    R = 8.31447;        % J/mol-K
    h = 6.62607e-34;    % J-s
    F = e * N_A;        % Coul/mole-of-unit-charge
    k_B = R / N_A;      % J/particle-K

    % Cantera gas object initialization
    gas = Solution('GRI30.yaml');
    H2 = speciesIndex(gas, 'H2');
    H2O = speciesIndex(gas, 'H2O');
    O2 = speciesIndex(gas, 'O2');
    N2 = speciesIndex(gas, 'N2');
    Nsp = nSpecies(gas);

    % Cell conditions
    Tcell = T_set;  % Temperature in K
    Pcell = P_set;            % Pressure in Pa

    % Anode composition and state
    xanode = zeros(Nsp, 1);
    xanode(H2O) = 0.03;
    xanode(H2) = 1 - xanode(H2O);
    set(gas, 'Temperature', Tcell, 'Pressure', Pcell, 'MoleFractions', xanode);
    muanode = chemPotentials(gas) / 1000;

    % Cathode composition and state
    xcathode = zeros(Nsp, 1);
    xcathode(O2) = 0.21;
    xcathode(N2) = 1 - xcathode(O2);
    set(gas, 'Temperature', Tcell, 'Pressure', Pcell, 'MoleFractions', xcathode);
    mucathode = chemPotentials(gas) / 1000;

    % Equilibrium chemical potentials
    muH2a_eq = muanode(H2);
    muH2Oa_eq = muanode(H2O);
    muO2c_eq = mucathode(O2);
    muEa_eq = 0;  % Electrochemical potential of an electron in the anode
    muOa_eq = muH2Oa_eq + 2 * muEa_eq - muH2a_eq;
    muOc_eq = muOa_eq;
    muEc_eq = (2 / 4) * muOc_eq - (1 / 4) * muO2c_eq;

    % Reaction kinetics and transport properties
    ioc = 1000;  % Cathode exchange current density in A/m^2
    Roc = ioc / (2 * F);
    Anode_Cathode_Ratio = 100;
    Roa = Anode_Cathode_Ratio * Roc;
    kappa_YSZ = 15;  % Ionic conductivity of YSZ in S/m
    kappa_O = kappa_YSZ / (-2 * F)^2;
    kappa_Ni = 2e6;  % Electronic conductivity of nickel in S/m
    kappa_E = kappa_Ni / (-1 * F)^2;

    % Geometric and transport parameters
    L_YSZ = 50e-6;  % Thickness of the YSZ electrolyte in meters
    L_GDLa = 5e-3;  % Thickness of the gas diffusion layer at anode in meters
    L_GDLc = 5e-3;  % Thickness of the gas diffusion layer at cathode in meters

    % Diffusivities
    concentration = Pcell / (R * Tcell);
    D_H2_H2O = Binary_Diffusivity_ijTP(iH2, iH2O, Tcell, Pcell);
    D_O2_N2 = Binary_Diffusivity_ijTP(iO2, iN2, Tcell, Pcell);

    % Iteration setup
    imax = 100 * 1000;  % Max current density in A/m^2
    vmax = imax / (2 * F);  % Max reaction velocity density in rxn/s/m^2
    vmin = 0;
    steps = 10000;
    dv = (vmax - vmin) / steps;
    v = vmin;
    found_voltage = false;
    ii = zeros(steps, 1);
    phiti = zeros(steps, 1);
    poweri = zeros(steps, 1);

    for idx = 1:steps
        % Iterative calculation code from your provided script
        JH2a = v;   % Flux, mol-i/m2-s
        Delta_xH2a = JH2a * L_GDLa / (concentration * D_H2_H2O);
        xH2a = xH2a_eq - Delta_xH2a;
        if (xH2a <= 0)
            disp('All H2 consumed...');
            break;
        end
        Delta_muH2a = -R * Tcell * log(1 - JH2a * L_GDLa / (xH2a_eq * concentration * D_H2_H2O));
        muH2a = muH2a_eq - Delta_muH2a;

        JH2Oa = v;  % Flux, mol-i/m2-s
        Delta_xH2Oa = JH2Oa * L_GDLa / (concentration * D_H2_H2O);
        xH2Oa = xH2Oa_eq + Delta_xH2Oa;
        if (xH2Oa >= 1)
            disp('Anode H2O hit unity...');
            break;
        end
        Delta_muH2Oa = R * Tcell * log(1 + JH2Oa * L_GDLa / (xH2Oa_eq * concentration * D_H2_H2O));
        muH2Oa = muH2Oa_eq + Delta_muH2Oa;

        JEa = 2 * v;  % Electron flux, mol-i/m2-s
        Delta_muEa = JEa * L_Nia / kappa_E;
        muEa = muEa_eq + Delta_muEa;

        muOa = muOa_eq + Delta_muH2a + R * Tcell * log(v / Roa + exp((Delta_muH2Oa + 2 * Delta_muEa) / (R * Tcell)));
        JOac = v;  % Oxygen ion flux, mol-i/m2-s
        Delta_muOac = JOac * L_YSZ / kappa_O;
        muOc = muOa + Delta_muOac;

        JO2c = 0.5 * v;  % Oxygen flux, mol-i/m2-s
        xO2c = 1 - (1 - xO2c_eq) * exp(JO2c * L_GDLc / (concentration * D_O2_N2));
        if (xO2c <= 0)
            disp('All O2 consumed...');
            break;
        end
        Delta_muO2c = -R * Tcell * log(xO2c / xO2c_eq);
        muO2c = muO2c_eq - Delta_muO2c;

        muEc = muEc_eq + 0.25 * Delta_muO2c + 0.5 * R * Tcell * log(v / Roc + exp((muOc - muOc_eq) / (R * Tcell)));

        JEc = 2 * v;  % Electron flux, mol-i/m2-s
        Delta_muEc = JEc * L_Nic / kappa_E;
        muEct = muEc + Delta_muEc;

        % Update the terminal voltage
        terminal_voltage = (muEct - muEat) / (-1 * F);

        % Save the calculated values
        ii(idx) = 2 * F * v;  
        phiti(idx) = terminal_voltage;
        poweri(idx) = ii(idx) * phiti(idx);

        % Check if the current terminal voltage matches the specified voltage
        if abs(terminal_voltage - specified_voltage) < 1e-2  % 0.01 V tolerance
            current_density = 2 * F * v;
            found_voltage = true;
            break;
        end

        v = v + dv;
    end

    if ~found_voltage
        disp('Specified voltage not reached within the current density range.');
        current_density = NaN;
        terminal_voltage = NaN;
        Pmax = NaN;
        return;
    end

    % Find the max power value and location
    [Pmax, ~] = max(poweri);  % Maximum power
end

