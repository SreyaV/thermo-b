function fugacity = fi_icrT(i,c,r,T)
% Return the fugacity (Pa) of species i in a mixture of composition 
% (molefractions) c, with density (kg/m3) r and temperature (K) T.
% C.F. Edwards, 1/11/09

% Use numerical differentiation and the fact that the fugacity is
% the partial of the normalized residual (including excess) Helmholtz
% function wrt moles of species i with fixed V, T, and other moles.

% Treat the composition vector as though it were mole numbers.  We will
% work with one kmol of the original mixture.
Mc = M_c(c);        % Molar mass.
Rc = R_c(c);        % Gas constant.
N = 1;              % Total mole number.
m = N*Mc;           % Mass of one kmol.
V = m/r;            % Volume of one kmol of mixture.
P = r*Rc*T;         % Ideal gas pressure for this density.

% Choose an increment in the mole number of species i.  Add/subtract it
% from the mixture.  Be careful of the boundaries.
Ninc = 1e-6;        % Try a ppm.

% Get the reduced temperature and density.
t = Tred_c(c)/T;
rred = rred_c(c);
d = r/rred;

% Start with the case where there is none of the species of interest.
if(c(i) == 0)                               % Must use forward difference.
    Narlow = N*ar_cdt(c,d,t);               % Norm. residual Helmholtz function.

    Nimid = c;                              % Copy the mole numbers.
    Nimid(i) = Ninc;                        % Add an increment of test species.
    Nmid   = N + Ninc;                      % Total mole number.
    cmid   = Nimid/Nmid;                    % New mole fractions.
    Mmid   = M_c(cmid);                     % Molar mass of mixture.
    mmid   = Mmid*Nmid;                     % Mass of mixture.
    rmid   = mmid/V;                        % Must hold volume fixed!
    rred   = rred_c(cmid);                  % Reducing value of density.
    dmid   = rmid/rred;                     % Reduced density.
    Narmid = Nmid*ar_cdt(cmid,dmid,t);      % N * ar product.
    
    Nihigh = c;                             % Copy the mole numbers.
    Nihigh(i) = 2*Ninc;                     % Add two increments of test species.
    Nhigh   = N + 2*Ninc;                   % Total mole number.   
    chigh   = Nihigh/Nhigh;                 % New mole fractions.
    Mhigh   = M_c(chigh);                   % Molar mass of mixture.
    mhigh   = Mhigh*Nhigh;                  % Mass of mixture.
    rhigh   = mhigh/V;                      % Must hold volume fixed!
    rred   = rred_c(chigh);                 % Reducing value of density.
    dhigh   = rhigh/rred;                   % Reduced density.
    Narhigh = Nhigh*ar_cdt(chigh,dhigh,t);  % N * ar product.

    % Use a second-order finite-difference formula.  See Ferziger's book on
    % numerical methods for details.
    dNardNi = (-3*Narlow + 4*Narmid - 1*Narhigh)/(2*Ninc); 
    fugacity = P*exp(dNardNi);
    return
end

% Now do the case where it is all the species of interest.
if(c(i) == 1)                               % Must use backward difference.
    Narhigh = N*ar_cdt(c,d,t);              % Norm. residual Helmholtz function.

    Nimid = c;                              % Copy the mole numbers.
    Nimid(i) = 1 - Ninc;                    % Remove an increment of test species.
    Nmid   = N - Ninc;                      % Total mole number.
    cmid   = Nimid/Nmid;                    % New mole fractions.
    Mmid   = M_c(cmid);                     % Molar mass of mixture.
    mmid   = Mmid*Nmid;                     % Mass of mixture.
    rmid   = mmid/V;                        % Must hold volume fixed!
    rred   = rred_c(cmid);                  % Reducing value of density.
    dmid   = rmid/rred;                     % Reduced density.
    Narmid = Nmid*ar_cdt(cmid,dmid,t);      % N * ar product.
    
    Nilow = c;                              % Copy the mole numbers.
    Nilow(i) = 1 - 2*Ninc;                  % Remove two increments of test species.
    Nlow   = N - 2*Ninc;                    % Total mole number.   
    clow   = Nilow/Nlow;                    % New mole fractions.
    Mlow   = M_c(clow);                     % Molar mass of mixture.
    mlow   = Mlow*Nlow;                     % Mass of mixture.
    rlow   = mlow/V;                        % Must hold volume fixed!
    rred   = rred_c(clow);                  % Reducing value of density.
    dlow   = rlow/rred;                     % Reduced density.
    Narlow = Nlow*ar_cdt(clow,dlow,t);      % N * ar product.

    % Use a second-order finite-difference formula.  See Ferziger's book on
    % numerical methods for details.
    dNardNi = (Narlow - 4*Narmid + 3*Narhigh)/(2*Ninc);
    fugacity = P*exp(dNardNi);
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
Nlow   = N - Ninc;                      % New total moles.
clow   = Nilow/Nlow;                    % New mole fractions.
Mlow   = M_c(clow);                     % Molar mass of mixture.
mlow   = Mlow*Nlow;                     % Mass of mixture.
rlow   = mlow/V;                        % Must hold volume fixed!
rred   = rred_c(clow);                  % Reducing value of density.
dlow   = rlow/rred;                     % Reduced density.
Narlow = Nlow*ar_cdt(clow,dlow,t);      % Helmholtz function.

Nihigh = c;                             % Same but high side.
Nihigh(i) = Nihigh(i) + Ninc;
Nhigh   = N + Ninc;
chigh   = Nihigh/Nhigh;
Mhigh   = M_c(chigh);
mhigh   = Mhigh*Nhigh;
rhigh   = mhigh/V;                         
rred    = rred_c(chigh);                % Reducing value of density.
dhigh   = rhigh/rred;
Narhigh = Nhigh*ar_cdt(chigh,dhigh,t);

% Use a second-order central difference formula.
dNardNi = (Narhigh - Narlow)/(2*Ninc);
fugacity = P*exp(dNardNi);
