% Mixture model for CO2/H2O of Paulus and Penoncello using the Lemmon
% and Jacobsen approach.
% Int. J. of Thermophysics, Vol. 27, p. 1373, 2006
% C.F. Edwards, 8/6/2007

% Values from published paper:

% Beta is one of the mole fraction exponents for temperature reducing.
Beta(CO2,H2O) = 0.974144949;
Beta(H2O,CO2) = Beta(CO2,H2O);

% Phi is the other mole fraction exponents for temperature reducing.
Phi(CO2,H2O) = 1;
Phi(H2O,CO2) = Phi(CO2,H2O);

% Zeta is the binary interaction coefficient for temperature reducing.
% It has units of Kelvin.
Zeta(CO2,H2O) = -302.4915666;
Zeta(H2O,CO2) = Zeta(CO2,H2O);

% Xi is the binary interaction coefficient for density reducing.
% It has units of molar density (kmols/m3)
Xi(CO2,H2O) = 2.4025883e-2;
Xi(H2O,CO2) = Xi(CO2,H2O);

% Fmix is the binary interaction coefficient for the excess Helmholtz
% energy.
Fmix(CO2,H2O) = 5.461671501;
Fmix(H2O,CO2) = Fmix(CO2,H2O);


% % Values from Steve Penoncello by email on 8/21/07:
% % Beta is the mole fraction exponent for temperature reducing.
% Beta(CO2,H2O) = 0.810570331;
% Beta(H2O,CO2) = Beta(CO2,H2O);
% % 
% % Phi is the other mole fraction exponents for temperature reducing.
% Phi(CO2,H2O) = 1;
% Phi(H2O,CO2) = Phi(CO2,H2O);
% 
% % Zeta is the binary interaction coefficient for temperature reducing.
% % It has units of Kelvin.
% Zeta(CO2,H2O) = -269.8733171;
% Zeta(H2O,CO2) = Zeta(CO2,H2O);
% 
% % Xi is the binary interaction coefficient for density reducing.
% % It has units of molar density (kmols/m3)
% Xi(CO2,H2O) = 1.4618624E-02;
% Xi(H2O,CO2) = Xi(CO2,H2O);
% 
% % Fmix is the binary interaction coefficient for the excess Helmholtz
% % energy.
% Fmix(CO2,H2O) = 5.755398892;
% Fmix(H2O,CO2) = Fmix(CO2,H2O);


% The general mixture model below is from Lemmon and Jacobsen
% Int. J. of Thermophysics, Vol. 20, p. 825, 1999

% Mixture excess fitting parameters:
Mix_i = [1 1 1 2 3 4 5 6 6 8];
Mix_j = [2 4 -2 1 4 4 4 0 4 -2];
Mix_N = [
    -0.245476271425e-1
    -0.241206117483
    -0.513801950309e-2
    -0.239824834123e-1
     0.259772344008
    -0.172014123104
     0.429490028551e-1
    -0.202108593862e-3
    -0.382984234857e-2
     0.262992331354e-5
    ];

% Note that the last Mix_N coefficient in the array above is different in
% the Paulus and Penoncello paper than in the Lemmon and Jacobsen paper.
% They give it as: 0.269923313540e-5 (a 2 is missing in the third place) as
% opposed to 0.262992331354e-5 in the original mixing model.
