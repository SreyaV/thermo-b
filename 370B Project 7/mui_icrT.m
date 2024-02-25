function dAdNi = mui_icrT(i,c,r,T)
% Return the chemical potential (J/kmol-species-i) of species i in a
% mixture of composition (molefractions) c, with density (kg/m3) r and
% temperature (K) T.
% C.F. Edwards, 2-16-10

% Use numerical differentiation and the fact that the chemical potential is
% the partial molar Helmholtz function with fixed V, T, and other moles.

% Treat the composition vector as though it were mole numbers.  We will
% work with one kmol of the original mixture.
N = 1;              % Total mole number.
m = N*M_c(c);       % Mass of one kmol.
V = m/r;            % Volume of one kmol of mixture.

% Choose an increment in the mole number of species i.  Add/subtract it
% from the mixture.  Be careful of the boundaries.
Ninc = 1e-6;        % Try a ppm.

% Start with the case where there is none of the species of interest.
if(c(i) == 0)                               % Must use forward difference.
    Alow = N*M_c(c)*a_crT(c,r,T);           % Helmholtz function.

    Nimid = c;                              % Copy the mole numbers.
    Nimid(i) = Ninc;                        % Add an increment of test species.
    Nmid  = N + Ninc;                       % Total mole number.
    cmid  = Nimid/Nmid;                     % New mole fractions.
    Mmid  = M_c(cmid);                      % Molar mass of mixture.
    mmid  = Mmid*Nmid;                      % Mass of mixture.
    rmid  = mmid/V;                         % Must hold volume fixed!
    Amid  = Nmid*Mmid*a_crT(cmid,rmid,T);
    
    Nihigh = c;                             % Copy the mole numbers.
    Nihigh(i) = 2*Ninc;                     % Add two increments of test species.
    Nhigh = N + 2*Ninc;                     % Total mole number.   
    chigh = Nihigh/Nhigh;                   % New mole fractions.
    Mhigh = M_c(chigh);                     % Molar mass of mixture.
    mhigh = Mhigh*Nhigh;                    % Mass of mixture.
    rhigh = mhigh/V;                        % Must hold volume fixed!
    Ahigh = Nhigh*Mhigh*a_crT(chigh,rhigh,T);

    % Use a second-order finite-difference formula.  See Ferziger's book on
    % numerical methods for details.
    dAdNi = (-3*Alow + 4*Amid - 1*Ahigh)/(2*Ninc); 
    return
end

% Now do the case where it is all the species of interest.
if(c(i) == 1)                               % Must use backward difference.
    Ahigh = N*M_c(c)*a_crT(c,r,T);          % Helmholtz function.

    Nimid = c;                              % Copy the mole numbers.
    Nimid(i) = 1 - Ninc;                    % Remove an increment of test species.
    Nmid  = N - Ninc;                       % Total mole number.
    cmid  = Nimid/Nmid;                     % New mole fractions.
    Mmid  = M_c(cmid);                      % Molar mass of mixture.
    mmid  = Mmid*Nmid;                      % Mass of mixture.
    rmid  = mmid/V;                         % Must hold volume fixed!
    Amid  = Nmid*Mmid*a_crT(cmid,rmid,T);
    
    Nilow = c;                              % Copy the mole numbers.
    Nilow(i) = 1 - 2*Ninc;                  % Remove two increments of test species.
    Nlow  = N - 2*Ninc;                     % Total mole number.   
    clow  = Nilow/Nlow;                     % New mole fractions.
    Mlow  = M_c(clow);                      % Molar mass of mixture.
    mlow  = Mlow*Nlow;                    % Mass of mixture.
    rlow  = mlow/V;                        % Must hold volume fixed!
    Alow  = Nlow*Mlow*a_crT(clow,rlow,T);

    % Use a second-order finite-difference formula.  See Ferziger's book on
    % numerical methods for details.
    dAdNi = (Alow - 4*Amid + 3*Ahigh)/(2*Ninc);
    return
end

% If asked for a value with a smaller amount of the species of interest,
% use the amount present as the increment.  (More restrictive.)
if(c(i) < Ninc)
    Ninc = c(i);
end

% If the increment would put us past unity molefraction for the species,
% scale it back.
if((c(i) + Ninc) > 1)
    Ninc = 1 - c(i);
end

% If we got to here we have a central difference case.
Nilow = c;                              % Copy the mole numbers.
Nilow(i) = Nilow(i) - Ninc;             % Decrement the test species.
Nlow = N - Ninc;                        % New total moles.
clow = Nilow/Nlow;                      % New mole fractions.
Mlow = M_c(clow);                       % Molar mass of mixture.
mlow = Mlow*Nlow;                       % Mass of mixture.
rlow = mlow/V;                          % Must hold volume fixed!
Alow = Nlow*Mlow*a_crT(clow,rlow,T);    % Helmholtz function.

Nihigh = c;                             % Same but high side.
Nihigh(i) = Nihigh(i) + Ninc;
Nhigh = N + Ninc;
chigh = Nihigh/Nhigh;
Mhigh = M_c(chigh);
mhigh = Mhigh*Nhigh;
rhigh = mhigh/V;                         
Ahigh = Nhigh*Mhigh*a_crT(chigh,rhigh,T);

% Use a second-order central difference formula.
dAdNi = (Ahigh - Alow)/(2*Ninc);
