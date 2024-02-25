function dmuvdci = dmuvdci_icTPrmu(i,c,T,P,rv,mu)
% Return the derivative of the chemical potential (J/kmol-species-i) of species i in a
% mixture of composition (molefractions) c, with density (kg/m3) r and
% temperature (K) T with respect to the molefraction of species i holding
% the moles of other species (and therefore the ratios of their
% molefractions) constant.
% C.F. Edwards, 2-16-10

global toler

% Use numerical differentiation and the fact that the chemical potential is
% the partial molar Helmholtz function with fixed V, T, and other moles.

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
    disp('Molefractions do not sum to unity in dmuvdci')
    return
end

% Watch out for pure substances.  You have no way to know which species
% should be REMOVED to form the derivative.
if(c(i) == 1)
    disp('Cannot form derivative for pure substance in dmuvdci')
    return
end

% Start with the case where there is little to none of the species of
% interest.
if(c(i) < cinc)                             % Use a forward difference.
    % Use mu passed from parent function for glow
    Nimid = c;                              % Copy the mole numbers.
    Ninc = cinc/(1-c(i)-cinc);              % Find the actual moles needed.
    Nimid(i) = Nimid(i) + Ninc;             % Increment the test species.
    Nmid = 1 + Ninc;                        % New total mole number.
    cmid = Nimid/Nmid;                      % New mole fractions.
    rv = rv_cTP(cmid,T,P,rv);
    gmid = mui_icrT(i,cmid,rv,T);           % Chemical potential function.

    Nihigh = c;                             % Copy the mole numbers.
    Ninc = cinc/(1-c(i)-2*cinc);            % Find the actual moles needed.
    Nihigh(i) = Nihigh(i) + Ninc;           % Increment the test species.
    Nhigh = 1 + Ninc;                       % New total mole number.
    chigh = Nihigh/Nhigh;                   % New mole fractions.
    rv = rv_cTP(chigh,T,P,rv);
    ghigh = mui_icrT(i,chigh,rv,T);         % Chemical potential function.

    % Use a second-order finite-difference formula.  See Ferziger's book on
    % numerical methods for details.
    dmuvdci = (4*gmid - ghigh - 3*mu)/(2*cinc);
    return
end

% Watch out for mixtures that are within cinc of being a pure substance.
% Adjust the increment down to accommodate.
if((c(i) + cinc) > 1)
    cinc = 1 - c(i);
end

% If we got to here we have a forward difference case.
% Use mu passed from parent function for g
Nihigh = c;                             % Copy the mole numbers.
Ninc = cinc/(1-c(i)-cinc);              % Find the actual moles needed.
Nihigh(i) = Nihigh(i) + Ninc;           % Increment the test species.
Nhigh = 1 + Ninc;                       % New total mole number.
chigh = Nihigh/Nhigh;                   % New mole fractions.
rv = rv_cTP(chigh,T,P,rv);
ghigh = mui_icrT(i,chigh,rv,T);         % Chemical potential function.

% Use a forward difference formula.
dmuvdci = (ghigh - mu)/(cinc);

