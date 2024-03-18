function Dij = Binary_Diffusivity_ijTP(i,j,T,P)
% Calculate ideal gas binary diffusivity for species pair i,j using 
% Eq. 16.3-1 from Bird, Stewart, and Lightfoot.
% C.F. Edwards, 2-29-08

Pa_atm = 101325;  % Pascals/atm

% Molar masses (kg/kmol), and critical props (K, Pa) for some species.
         N2 = 1;  O2 = 2;  Ar = 3; CO2 = 4;  H2O = 5;  CO = 6;  H2 = 7; CH4 = 8;
M_i  = [ 28.014;  31.999;  39.948;  44.010;   18.015;  28.010;   2.016;  16.043];
Tc_i = [  126.2;   154.6;   150.9;   304.2;    647.1;   132.9;   33.19;   190.6];
Pc_i = [34.00e5; 50.43e5; 48.98e5; 73.76e5; 220.55e5; 34.99e5; 13.13e5; 45.99e5];

% Change units on the pressure.
P_atm = P/Pa_atm;

% Set the model constants.
with_water = (i == 5)||(j == 5);
if(with_water)
    a = 3.640e-4;
    b = 2.334;
else
    a = 2.745e-4;
    b = 1.823;
end

% Set the molar masses.
Mi = M_i(i);
Mj = M_i(j);

% Set the critical pressures.
Pc_atm_i = Pc_i(i)/Pa_atm;
Pc_atm_j = Pc_i(j)/Pa_atm;

% Set the critical temperatures.
Tcrit_i = Tc_i(i);
Tcrit_j = Tc_i(j);

% Note that Eq. 16.3-1 has imposed units as follows:
% P -> atm
% T -> K
% M -> kg/kmol
% Dij -> cm2/s
% Get the pressure-diffusivity product from Eq. 16.3-1
PDij = (Pc_atm_i*Pc_atm_j)^(1/3) *...
    (Tcrit_i*Tcrit_j)^(5/12) *...
    (1/Mi + 1/Mj)^(1/2) *...
    a*(T/sqrt(Tcrit_i*Tcrit_j))^b;

% Get the diffusivity (in cm2/s).
Dij = PDij/P_atm;

% Change the units to m2/s.
Dij = Dij/(100^2);
return

