function dfvdci = dfvdci_icTPrf(i,c,T,P,rv,f)
% Return the derivative of the fugacity (Pa) of species i in a
% mixture of composition (molefractions) c, with density (kg/m3) r and
% temperature (K) T with respect to the molefraction of species i holding
% the moles of other species (and therefore the ratios of their
% molefractions) constant.
% C.F. Edwards, 2-16-10

global toler

% Use numerical differentiation and the fact that the fugacity is
% the partial derivative of the product of mole number and normalized
% Helmholtz function with species i mole number at fixed V, T, and other moles.

% Treat the composition vector as though it were mole numbers.  We will
% work with one kmol of the original mixture.

% Choose an increment in the molefraction of species i.  Add/subtract the
% required increment of moles for this amount of molefraction change to and
% from the mixture.  Be careful of the boundaries.
cinc = 1e-6;        % Try a ppm.

% Check to make sure you have a viable composition.
csum = sum(c);
if(abs(csum-1) > toler)
    c
    disp('Molefractions do not sum to unity in dfvdci')
    return
end

% Watch out for pure substances.  You have no way to know which species
% should be REMOVED to form the derivative.
if(c(i) == 1)
    disp('Cannot form derivative for pure substance in dfvdci')
    return
end

% Start with the case where there is little to none of the species of
% interest.
if(c(i) < cinc)                             % Use a forward difference.
    % Use f passed from parent function for flow
    Nimid = c;                              % Copy the mole numbers.
    Ninc = cinc/(1-c(i)-cinc);              % Find the actual moles needed.
    Nimid(i) = Nimid(i) + Ninc;             % Increment the test species.
    Nmid = 1 + Ninc;                        % New total mole number.
    cmid = Nimid/Nmid;                      % New mole fractions.
    rv = rv_cTP(cmid,T,P,rv);               % Vapor density.
    fmid = fi_icrT(i,cmid,rv,T);            % Fugacity.

    Nihigh = c;                             % Copy the mole numbers.
    Ninc = cinc/(1-c(i)-2*cinc);            % Find the actual moles needed.
    Nihigh(i) = Nihigh(i) + Ninc;           % Increment the test species.
    Nhigh = 1 + Ninc;                       % New total mole number.
    chigh = Nihigh/Nhigh;                   % New mole fractions.
    rv = rv_cTP(chigh,T,P,rv);              % Vapor density.
    fhigh = fi_icrT(i,chigh,rv,T);          % Fugacity.

    % Use a second-order finite-difference formula.  See Ferziger's book on
    % numerical methods for details.
    dfvdci = (4*fmid - fhigh - 3*f)/(2*cinc);
    return
end

% Watch out for mixtures that are within cinc of being a pure substance.
% Adjust the increment down to accommodate.
if((c(i) + cinc) > 1)
    cinc = 1 - c(i);
end

% If we got to here we have a forward difference case.
% Use f passed from parent function as flow.
Nihigh = c;                             % Copy the mole numbers.
Ninc = cinc/(1-c(i)-cinc);              % Find the actual moles needed.
Nihigh(i) = Nihigh(i) + Ninc;           % Increment the test species.
Nhigh = 1 + Ninc;                       % New total mole number.
chigh = Nihigh/Nhigh;                   % New mole fractions.
rv = rv_cTP(chigh,T,P,rv);              % Vapor density.
fhigh = fi_icrT(i,chigh,rv,T);          % Fugacity.

% Use a forward difference formula.
dfvdci = (fhigh - f)/(cinc);
