function [i mu xac] = SOFC_Element_V(V,x_eq,mu_eq,Tcell,K,L,ioa,ioc)
% Return the current density i (A/m2) and the array of electrochemical 
% potentials (J/mol) for a differential area element of a hydrogen-water-air 
% SOFC operating at electrical potential difference between the cathode 
% and anode (V) with channel (and therefore equilibrium) 
% chemical potentials in array mu_eq (J/mol), conductivities K 
% (mol-i/m2/s)/(J/mol-i), path lengths L (m), and anode and cathode
% exchange current densities ioa and ioc (A/m2/s).

% Newton method

eps = 1e-4;
delta = 1;
iterations = 0;

% initial guess
i0 = 1; %A/m²

while delta>eps
    
    iterations = iterations +1;

    % Compute diff
    di = 0.1; %A/m²
    [V_k mu_k xac_k] = SOFC_Element_icTKL(i0,x_eq,mu_eq,Tcell,K,L,ioa,ioc);
    [V_k_dk mu_k_dk xac_k_dk] = SOFC_Element_icTKL(i0+di,x_eq,mu_eq,Tcell,K,L,ioa,ioc);
    diff = (V_k_dk-V_k) /di;

    % Update
    i0 = i0 - V_k/diff;
    delta = abs(V-V_k);

    if i0<0
        i0 = 0;
    end

end

% Output
i = i0;
mu = mu_k;
xac = xac_k;