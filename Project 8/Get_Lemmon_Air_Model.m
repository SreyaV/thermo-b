% Mixture model for N2/O2/Ar of Lemmon et al.
% J. Phys. Chem. Ref. Data, Vol. 29, p. 331, 2000

Ru = 8314.510;  % J/kmol_K

% Beta is one of the mole fraction exponents for temperature reducing.
Beta(N2,O2) = 1;
Beta(O2,N2) = Beta(N2,O2);
Beta(N2,Ar) = 1;
Beta(Ar,N2) = Beta(N2,Ar);
Beta(O2,Ar) = 1;
Beta(Ar,O2) = Beta(O2,Ar);

% Phi is the other mole fraction exponents for temperature reducing.
Phi(N2,O2) = 1;
Phi(O2,N2) = Phi(N2,O2);
Phi(N2,Ar) = 1;
Phi(Ar,N2) = Phi(N2,Ar);
Phi(O2,Ar) = 1;
Phi(Ar,O2) = Phi(O2,Ar);

% Zeta is the binary interaction coefficient for temperature reducing.
Zeta(N2,O2) = -0.856350;
Zeta(O2,N2) = Zeta(N2,O2);
Zeta(N2,Ar) = -1.237713;
Zeta(Ar,N2) = Zeta(N2,Ar);
Zeta(O2,Ar) = -2.115126;
Zeta(Ar,O2) = Zeta(O2,Ar);

% Xi is the binary interaction coefficient for density reducing.
Xi(N2,O2) = -0.00041847;
Xi(O2,N2) = Xi(N2,O2);
Xi(N2,Ar) = -0.00076031;
Xi(Ar,N2) = Xi(N2,Ar);
Xi(O2,Ar) =  0.00041232;
Xi(Ar,O2) = Xi(O2,Ar);

% Fmix is the binary interaction coefficient for the excess Helmholtz
% energy.
Fmix(N2,O2) = 1;
Fmix(O2,N2) = Fmix(N2,O2);
Fmix(N2,Ar) = 1.121527;
Fmix(Ar,N2) = Fmix(N2,Ar);
Fmix(O2,Ar) = 0.597203;
Fmix(Ar,O2) = Fmix(O2,Ar);

% Mixture excess fitting parameters:
Mix_i = [2 2];
Mix_j = [-1.4 1.5];
Mix_N = [-0.00195245 0.00871334];
