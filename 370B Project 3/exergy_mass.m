function x = exergy_mass(fluid)
% Return the mass-specific internal exergy of the fluid (J/kg);
% C.F. Edwards, 1-10-10

global To Po mu_o
global g_tm

% The variable g_tm is used to save the Gibbs function of a pure fluid in
% the thermomechanical dead state.  This saves time recalculating it.  Note
% that it will not be correct if you change fluids.  So if you do, reset
% its value to zero in the main routine to force its recalculation.

% If g_tm did not get initialized outside of this routine, do so by  
% forcing its calculation.
if isempty(g_tm)
    g_tm = 0;
end

% Get the basics from the fluid.
u = intEnergy_mass(fluid);
r = density(fluid);
s = entropy_mass(fluid);
a = u + Po/r - To*s;

% Figure out which fluid you have and deal with it accordingly.
% Use the first/only species name to check.
species1 = char(speciesName(fluid,1));

% Check two character names.
if((nSpecies(fluid) == 1) && (length(species1) == 2))
    if(strcmp(species1,'H2'))
        if(g_tm == 0)   % Don't calc g_tm if it is already set
            ref = Solution('liquidvapor.xml','hydrogen');
            set(ref,'T',To,'P',Po);
            g_tm = gibbs_mass(ref);
        end
        xchem = 1.165e+008;
        x = a - g_tm + xchem;
        return
    end
    
    if(strcmp(species1,'N2'))
        if(g_tm == 0)   % Don't calc g_tm if it is already set
            ref = Solution('liquidvapor.xml','nitrogen');
            set(ref,'T',To,'P',Po);
            g_tm = gibbs_mass(ref);
        end
        xchem = 2.503e+004;
        x = a - g_tm + xchem;
        return
    end
    
    if(strcmp(species1,'O2'))
        if(g_tm == 0)   % Don't calc g_tm if it is already set
            ref = Solution('liquidvapor.xml','oxygen');
            set(ref,'T',To,'P',Po);
            g_tm = gibbs_mass(ref);
        end
        xchem = 1.239e+005;
        x = a - g_tm + xchem;
        return
    end
end

% Check three character names.
if((nSpecies(fluid) == 1) && (length(species1) == 3))
    if(strcmp(species1,'H2O'))
        if(g_tm == 0)   % Don't calc g_tm if it is already set
            ref = Solution('liquidvapor.xml','water');
            set(ref,'T',To,'P',Po);
            g_tm = gibbs_mass(ref);
        end
        % Assume the dead state is saturated so xchem = 0.
        x = a - g_tm;
        return
    end
    
    if(strcmp(species1,'CO2'))
        if(g_tm == 0)   % Don't calc g_tm if it is already set
            ref = Solution('liquidvapor.xml','carbondioxide');
            set(ref,'T',To,'P',Po);
            g_tm = gibbs_mass(ref);
        end
        xchem = 4.456e+005;
        x = a - g_tm + xchem;
        return
    end
end

% Check six character names.
if((nSpecies(fluid) == 1) && (length(species1) == 6))
    if(strcmp(species1,'C2F4H2'))
        if(g_tm == 0)   % Don't calc g_tm if it is already set
            ref = Solution('liquidvapor.xml','hfc134a');
            set(ref,'T',To,'P',Po);
            g_tm = gibbs_mass(ref);
        end
        % R134a does not react to environmental species, so set xchem = 0.
        x = a - g_tm;
        return
    end
    % If we ever have other 6 char names, put them here!
end

% If you get down to here, it is the GRI30 property set so do full chemical
% exergy.  Note that the dead state chemical potentials must already be
% loaded into the global array mu_o.  (If you calculated them here,
% repeated calls would take forever.)

% When a non-environmental species is converted it does so by the reaction:
% CxHyOzNwARq -> xCO2 + (y/2)H2O + (w/2)N2 + qAR - (x + y/4 - z/2)O2
% Note that this also treats environmental species present in the resource
% correctly--they are converted to themselves by the reaction!

N  = moleFractions(fluid);
M  = meanMolecularWeight(fluid);
Ns = nSpecies(fluid);

% Define indices to find the environmental species in the array.
% Note that argon as a species goes by "AR" while argon as an element goes
% by "Ar".  Isn't Cantera wonderful some times?
iCO2 = speciesIndex(fluid,'CO2');
iH2O = speciesIndex(fluid,'H2O');
iN2  = speciesIndex(fluid,'N2');
iO2  = speciesIndex(fluid,'O2');
iAR  = speciesIndex(fluid,'AR');

G_o = 0;
for i=1:1:Ns
    if(N(i) ~= 0.0)
        % function n = nAtoms(a,k,m)
        % NATOMS-Number of atoms of m in species k.
        N_O  = nAtoms(fluid,i,elementIndex(fluid,'O'));
        N_H  = nAtoms(fluid,i,elementIndex(fluid,'H'));
        N_C  = nAtoms(fluid,i,elementIndex(fluid,'C'));
        N_N  = nAtoms(fluid,i,elementIndex(fluid,'N'));
        N_AR = nAtoms(fluid,i,elementIndex(fluid,'Ar'));
        
        N_CO2 = N_C;
        N_H2O = N_H/2;
        N_N2  = N_N/2;
        N_O2  = N_O/2 - N_C - N_H/4;
        
        G_o = G_o + N(i)*(N_CO2*mu_o(iCO2) + (N_H2O)*mu_o(iH2O) + (N_N2)*mu_o(iN2)...
            + (N_AR)*mu_o(iAR) + (N_O2)*mu_o(iO2));
    end
end
x = a - G_o/M;
