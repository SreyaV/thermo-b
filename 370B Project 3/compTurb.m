function [Pret,Tret,finalState,pathStates] = compTurb(fluid,P1,T1,P2,eta_p,nsteps)
% generalized function for a general fluid compressor or turbine
% performing polytropic compression/expansion
% Pret is the process path pressure that the compressor/turbine takes
% Tret is the process path temperature that the compressor/turbine takes
%
% [Pret,Tret,finalState,pathStates] = compTurb(fluid,P1,T1,P2,eta_p,nsteps)

% find out the pressure step size
dP = (P2 - P1)/nsteps;
Pret = P1:dP:P2;
Tret = 0*Pret;
hret = 0*Pret;
sret = 0*Pret;
vF = 0*Pret; 

% maybe don't do this if the fluid is supposed to be supplied at its
% initial condition directly into the function
try
    set(fluid,'T',T1,'P',P1);
catch % if the fluid is in its saturated state, set using its vapor fraction
    set(fluid,'P',P1,'Vapor',vaporFraction(fluid));
end

% set initial enthalpy and entropy
h1 = enthalpy_mass(fluid);
s1 = entropy_mass(fluid);

% loop over pressure
for iter = 1:1:length(Pret)
    % take the gas to the hypothetical isentropic state dP away from 1
    set(fluid,'S',s1,'P',Pret(iter));
    h2_s = enthalpy_mass(fluid);
    % use the definition of polytropic efficiency to get the actual h2
    if P2>=P1 % if compressing
        h2 = h1 + (h2_s-h1)/eta_p;
    else % if expanding
        h2 = h1 + (h2_s-h1)*eta_p;
    end
    set(fluid,'H',h2,'P',Pret(iter));
    % extract all required parameters
    Tret(iter) = temperature(fluid);
    hret(iter) = enthalpy_mass(fluid);
    sret(iter) = entropy_mass(fluid);
    % if the fluid is a fluid, get its vapor fraction
    try
        vF(iter) = vaporFraction(fluid);
    catch
        vF(iter) = nan;
    end
    h1 = enthalpy_mass(fluid);
    s1 = entropy_mass(fluid);
end

finalState.T = Tret(end);
finalState.P = Pret(end);
finalState.H = enthalpy_mass(fluid);
finalState.S = entropy_mass(fluid);

pathStates.T = Tret;
pathStates.P = Pret;
pathStates.H = hret;
pathStates.S = sret;
pathStates.vF = vF;

% find out exergy if required else forget about it
pathStates.X = 0;
pathStates.Xf = 0;


end

