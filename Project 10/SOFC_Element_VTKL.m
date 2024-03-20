function [ic,mu,xac,delta, phi] = SOFC_Element_VTKL(V,x_eq,mu_eq,Tcell,K,L,ioa,ioc,ic_guess,varagin)

% Return the current density and the array of electrochemical potentials (J/mol) for a
% differential area element of a hydrogen-water-air SOFC operating at
% current density i (A/m2) with channel (and therefore equilibrium) 
% chemical potentials in array mu_eq (J/mol), conductivities K 
% (mol-i/m2/s)/(J/mol-i), path lengths L (m), and anode and cathode
% exchange current densities ioa and ioc (A/m2/s).
% delta is the error in matching the voltage

global R F

if nargin == 1
    phi_eq = varagin{1};
    if V>phi_eq
        fprintf('\nYou are stoopid. Voltage asked is much too high!!!\n');
        ic = nan; mu = nan; xac = nan; delta = nan;
        return;
    end
end

helperFun = @(icIter)(SOFC_Element_icTKL(icIter,x_eq,mu_eq,Tcell,K,L,ioa,ioc)-V);

% options = optimoptions('fsolve','OptimalityTolerance',1e-12,'Display','none');
options = optimset('Display','none','TolFun',1e-12);
icSol = fzero(helperFun,ic_guess,options);

[phi,mu,xac] = SOFC_Element_icTKL(icSol,x_eq,mu_eq,Tcell,K,L,ioa,ioc);

ic = icSol;

delta = phi - V;

end