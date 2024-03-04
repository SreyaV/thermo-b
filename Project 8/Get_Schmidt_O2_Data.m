% Loads data into an array for dealing with multiple species.
% Data file for oxygen properties using the data from 
% "A New Form of the Equation of State for Pure Substances and its
% Application to Oxygen," Fluid Phase Equilibria, 19, 175 1985.
% C.F. Edwards, 9/2/07 

% Constants and ideal-gas ref. state properties:
Ru         = 8314.34;                       % J/kmol-K
M_i(O2)    = 31.9988;                       % kg/kmol
R_i(O2)    = Ru/M_i(O2);                    % J/kg-K
wP_i(O2)   = 0.022;                         % Pitzer's acentric factor
Tref_i(O2) = 298.15;                        % K
Pref_i(O2) = 101325;                        % Pa
rref_i(O2) = Pref_i(O2)/R_i(O2)/Tref_i(O2); % kg/m3

% Fixed-point properties and limits:
Tcrit_i(O2) = 154.581;                      % K
Pcrit_i(O2) = 5.043e6;                      % Pa
rcrit_i(O2) = 13.63*M_i(O2);                % kg/m3

Ttrip_i(O2)  = 54.361;                      % K
Ptrip_i(O2)  = 146.33;                      % Pa
rftrip_i(O2) = 40.816*M_i(O2);              % kg/m3
rgtrip_i(O2) = 0.0003237*M_i(O2);           % kg/m3

Tupper_i(O2) = 300;                         % K
Tlower_i(O2) = 54;                          % K
Pupper_i(O2) = 818e5;                       % Pa
rupper_i(O2) = 41*M_i(O2);                  % kg/m3

% Ideal gas Helmholtz fit coefficients for Eq.15.
FR_Npoly0(O2) = 4;
FR_Neinst(O2) = 0;
FR_Nspec0(O2) = 2;
FR_t0(1:4,O2) = [0; 1.5; -2; 1];
% The original Schmidt and Wagner paper uses zero values for h0ref and s0ref.
% The paper by Stewart et al. that gives tables from the Schmidt and Wagner
% fits however uses the reference values below:
href_i(O2) = 8682000./M_i(O2);              % J/kg
sref_i(O2) = 205037./M_i(O2);               % J/kg-K
% The Lemmon air model uses these values:
% href_i(O2) = 8680000./M_i(O2);              % J/kg
% sref_i(O2) = 205043./M_i(O2);               % J/kg-K
% Note that the order of the coefficients is different from those in the 
% original paper:  We put the log term after the polynomial terms and
% there are two special terms.  Here are the values from Table 3:
k1 = -0.740775e-3;
k2 = -0.664930e-4;
k3 =  0.250042e+1;
k4 = -0.214487e+2;
k5 =  0.101258e+1;
k6 = -0.944365;
k7 =  0.145066e+2;
k8 =  0.749148e+2;
k9 =  0.414817e+1;
d0 = rref_i(O2)/rcrit_i(O2);
t0 = Tcrit_i(O2)/Tref_i(O2);
FR_N0(1,O2) = k9 - sref_i(O2)/R_i(O2) - log(d0);
FR_N0(2:8,O2) = [k1; k2; k4 + href_i(O2)/(R_i(O2)*Tref_i(O2)*t0); k3; 0; k5; k6];
% The special-term coefficients as listed in Eq. 15 in the original paper.
FR_spec0(7:8,O2) = [k7; k8];

% Number of terms of each type in the residual part of the FR:
FR_Npoly(O2)  = 13; 
FR_Nexp(O2)   = 32-13;

% Fundamental Relation coefficients from Eq. 11 and Table 2 of the paper:
Table_2 = [
1    0.3983768749       1    0.0     0
2   -0.1846157454e+1    1    1.5     0
3    0.4183473197       1    2.5     0
4    0.2370620711e-1    2   -0.5     0
5    0.9771730573e-1    2    1.5     0
6    0.3017891294e-1    2    2.0     0
7    0.2273353212e-1    3    0.0     0
8    0.1357254086e-1    3    1.0     0
9   -0.4052698943e-1    3    2.5     0
10   0.5454628515e-3    6    0.0     0
11   0.5113182277e-3    7    2.0     0
12   0.2953466883e-6    7    5.0     0
13  -0.8687645072e-4    8    2.0     0
14  -0.2127082589       1    5.0     2
15   0.8735941958e-1    1    6.0     2
16   0.1275509190       2    3.5     2
17  -0.9067701064e-1    2    5.5     2
18  -0.3540084206e-1    3    3.0     2
19  -0.3623278059e-1    3    7.0     2
20   0.1327699290e-1    5    6.0     2
21  -0.3254111865e-3    6    8.5     2
22  -0.8313582932e-2    7    4.0     2
23   0.2124570559e-2    8    6.5     2
24  -0.8325206232e-3   10    5.5     2
25  -0.2626173276e-4    2   22.0     4
26   0.2599581482e-2    3   11.0     4
27   0.9984649663e-2    3   18.0     4
28   0.2199923153e-2    4   11.0     4
29  -0.2591350486e-1    4   23.0     4
30  -0.1259630848       5   17.0     4
31   0.1478355637       5   18.0     4
32  -0.1011251078e-1    5   23.0     4
 ];

FR_N(1:32,O2) = Table_2(:,2);
FR_d(1:32,O2) = Table_2(:,3);
FR_t(1:32,O2) = Table_2(:,4);
FR_c(1:32,O2) = Table_2(:,5);
