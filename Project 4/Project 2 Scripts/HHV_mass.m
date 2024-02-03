function HHV = HHV_mass(gas)
% Find the higher heating value of a gas.  Use extra oxygen to ensure
% complete combustion.  Since To is so low, dissociation is not an issue.
% Note that Cantera can have trouble with 25C so use 300K.  Ug!
% C.F. Edwards, 1-10-10

% Set a fractional tolerance value.
toler = 1e-6;

M    = molecularWeights(gas);
Nsp  = nSpecies(gas);
Nel  = nElements(gas);
iO2  = speciesIndex(gas,'O2');
iH2O = speciesIndex(gas,'H2O');
iCO2  = speciesIndex(gas,'CO2');
iSO2  = speciesIndex(gas,'SO2');

% Save the state so we can reset it at the end.
Tin = temperature(gas);
Pin = pressure(gas);
xin = moleFractions(gas);

% Set the temperature and pressure at which to evaluate the heating value.
% Sometimes Cantera will crash Matlab if you use 25C.  So you might want to
% use 300K instead.  It makes no difference to the heating value.
% Tref = 25+273.15;   % K
Tref = 300;         % K
Pref = 101325;      % Pa
% Use one mole of fuel mixture.
Nfuel = moleFractions(gas); 

% Cantera can have trouble equilibrating with small, nonzero values for
% trace species.  It will make no difference to the heating value if we
% kill them.
for j=1:1:Nsp
    if(Nfuel(j) < toler)
        Nfuel(j) = 0;
    end
end
% Renormalize and get the molecular mass.
Nfuel = Nfuel/sum(Nfuel);
set(gas,'T',Tref,'P',Pref,'X',Nfuel);
Mfuel = meanMolecularWeight(gas);

% Find out which elements the gas contains.  We only care about C,H,O,S.
Has_H = 0;
Has_C = 0;
Has_S = 0;
Has_O = 0;
for i=1:1:Nel
    if(strcmp(elementName(gas,i),'H'))
        Has_H = 1;
    elseif(strcmp(elementName(gas,i),'C'))
        Has_C = 1;
    elseif(strcmp(elementName(gas,i),'S'))
        Has_S = 1;
    elseif(strcmp(elementName(gas,i),'O'))
        Has_O = 1;
    end
end

% Find the stoichiometric oxygen requirment for this fuel.
N_H = 0;
N_C = 0;
N_S = 0;
N_O = 0;
for i=1:1:Nsp
    This_Species_Counts = (Nfuel(i)~=0.0)&&(i~=iH2O)&&(i~=iCO2)&&(i~=iSO2);
    if This_Species_Counts
        % function n = nAtoms(gas,i,m)
        % NATOMS-Number of atoms of m in species k.
        if(Has_H)
            m = elementIndex(gas,'H');
            N_H  = N_H + Nfuel(i)*nAtoms(gas,i,m);
        end
        if(Has_C)
            m = elementIndex(gas,'C');
            N_C  = N_C + Nfuel(i)*nAtoms(gas,i,m);
        end
        if(Has_S)
            m = elementIndex(gas,'S');
            N_S  = N_S + Nfuel(i)*nAtoms(gas,i,m);
        end
        if(Has_O)
            m = elementIndex(gas,'O');
            N_O  = N_O + Nfuel(i)*nAtoms(gas,i,m);
        end
    end
end
% We have the total number of each of the relevant atoms in the fuel.  Now
% find the oxygen needed by the H, C, and S atoms.  Note that if the fuel
% already contains some oxygen that will be present above and beyond the
% amount we choose to add.
O_needed   = N_H/2 + 2*N_C + 2*N_S;
O2_needed  = O_needed/2;

% Check to see if this is not combustible.
if O_needed == 0
    % Restore state of gas.
    set(gas,'T',Tin,'P',Pin,'X',xin);
    % Set the LHV to zero.
    HHV = 0;
    return
end

% Add enough oxygen to get 100% excess air.  This will be enough to
% suppress CO and H2.
Noxid = 2*O2_needed;
% Make a mixture.
Nmix = Nfuel;                   % Put fuel in the mixture
Nmix(iO2) = Nmix(iO2) + Noxid;  % Add oxygen in excess of existing

% Find the mass of the mix (fuel plus added oxid).
mass_mix  = 0;
for j=1:1:Nsp
    mass_mix  = mass_mix  + Nmix(j)*M(j);
end
% Since there is one mole of fuel in the mix, the molecular mass gives
% its contribution to the total mixture mass.
mass_fraction_fuel  = Mfuel/mass_mix;
% Normalize the mix.
Xmix = Nmix/sum(Nmix);
% Do the heating value problem.
set(gas,'T',Tref,'P',Pref,'X',Xmix);        % Put in the mixture
h_reactants = enthalpy_mass(gas);           % Find its enthalpy
equilibrate(gas,'TP');                      % Burn it at constant TP
h_products = enthalpy_mass(gas);            % Find the enthalpy
% Find the mass fraction of water.
y_products = massFractions(gas);
yH2O = y_products(iH2O);
% Find hfg for water.
water = Solution('liquidvapor.yaml','water');
setState_Tsat(water,[Tref 1]);
hg = enthalpy_mass(water);
setState_Tsat(water,[Tref 0]);
hf = enthalpy_mass(water);
hfgH2O = hg-hf;

% The enthalpy difference per unit mass of fuel is the heating value.
% Since the water is all liquid, this is the higher heating value.
HHV = (h_reactants - (h_products-yH2O*hfgH2O))/mass_fraction_fuel;

% Restore state of gas.
set(gas,'T',Tin,'P',Pin,'X',xin);
