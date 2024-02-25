% Loads data into an array for dealing with multiple species.
% Data file for argon properties using the data from 
% "A New Equation of State for Argon..." by Tegeler, Span, and Wagner,
% J. Phy. Chem. Ref. Data, 28(3), 779, 1999
% C.F. Edwards, 8/31/07 

% Constants and ideal-gas ref. state properties:
Ru         = 8314.51;                       % J/kmol-K
M_i(Ar)    = 39.948;                        % kg/kmol
R_i(Ar)    = 208.1333;                      % J/kg-K
wP_i(Ar)   = -0.00219;                      % Pitzer's acentric factor
Tref_i(Ar) = 298.15;                        % K
Pref_i(Ar) = 101325;                        % Pa
rref_i(Ar) = Pref_i(Ar)/R_i(Ar)/Tref_i(Ar); % kg/m3

% Fixed-point properties and limits:
Tcrit_i(Ar) = 150.687;                      % K
Pcrit_i(Ar) = 4.863e6;                      % Pa
rcrit_i(Ar) = 535.6;                        % kg/m3

Ttrip_i(Ar)  = 83.8058;                     % K
Ptrip_i(Ar)  = 68891;                       % Pa
rftrip_i(Ar) = 1416.77;                     % kg/m3
rgtrip_i(Ar) = 4.0546;                      % kg/m3

Tupper_i(Ar) = 700;                         % K
Tlower_i(Ar) = 83.8058;                     % K
Pupper_i(Ar) = 1000e6;                      % Pa
rupper_i(Ar) = 45.8136*M_i(Ar);             % kg/m3

% Ideal gas Helmholtz fit coefficients for Eq.4.6.
FR_Npoly0(Ar) = 2;
FR_Neinst(Ar) = 0;
FR_t0(1:2,Ar) = [0; 1];
% The original Tegeler et al. paper uses zero for h0ref and s0FRref.  Use these
% values and the associated FR_N0 values if you want to compare with the
% tables in the original paper.
% href_i(Ar) = 0;                           % J/kg
% sref_i(Ar) = 0;                           % J/kg-K
% FR_N0(1:3,Ar) = [8.31666243; -4.94651164; 1.5];
% Lemmon's air model uses the nonzero values below.  The value of FR_N0(2)
% determines the reference enthalpy while the value of FR_N0(1) then sets
% the entropy.  Be sure to adjust the enthalpy first since this constant
% affects the entropy too.
href_i(Ar) = 6197000./M_i(Ar);              % J/kg
sref_i(Ar) = 154737./M_i(Ar);               % J/kg-K
FR_N0(1:3,Ar) = [-10.2938169971; -3.40969594026e-4; 1.5];

% Number of terms of each type in the residual part of the FR:
FR_Npoly(Ar)  = 12; 
FR_Nexp(Ar)   = 37-12;
FR_Ngaus(Ar)  = 41-37;

% Fundamental Relation coefficients from Table 30 of the paper:
Table_30 = [
1    0.88722304990011e-1    1    0.00    0    0    0    0    0
2    0.70514805167298       1    0.25    0    0    0    0    0
3   -0.16820115654090e+1    1    1.00    0    0    0    0    0
4   -0.14909014431486       1    2.75    0    0    0    0    0
5   -0.12024804600940       1    4.00    0    0    0    0    0
6   -0.12164978798599       2    0.00    0    0    0    0    0
7    0.40035933626752       2    0.25    0    0    0    0    0
8   -0.27136062699129       2    0.75    0    0    0    0    0
9    0.24211924579645       2    2.75    0    0    0    0    0
10   0.57889583185570e-2    3    0.00    0    0    0    0    0
11  -0.41097335615341e-1    3    2.00    0    0    0    0    0
12   0.24710761541614e-1    4    0.75    0    0    0    0    0
13  -0.32181391750702       1    3.00    1    0    0    0    0
14   0.33230017695794       1    3.50    1    0    0    0    0
15   0.31019986287345e-1    3    1.00    1    0    0    0    0
16  -0.30777086002437e-1    4    2.00    1    0    0    0    0
17   0.93891137419581e-1    4    4.00    1    0    0    0    0
18  -0.90643210682031e-1    5    3.00    1    0    0    0    0
19  -0.45778349276654e-3    7    0.00    1    0    0    0    0
20  -0.82659729025197e-4   10    0.50    1    0    0    0    0
21   0.13013415603147e-3   10    1.00    1    0    0    0    0
22  -0.11397840001996e-1    2    1.00    2    0    0    0    0
23  -0.24455169960535e-1    2    7.00    2    0    0    0    0
24  -0.64324067175955e-1    4    5.00    2    0    0    0    0
25   0.58889471093674e-1    4    6.00    2    0    0    0    0
26  -0.64933552112965e-3    8    6.00    2    0    0    0    0
27  -0.13889862158435e-1    3   10.00    3    0    0    0    0
28   0.40489839296910       5   13.00    3    0    0    0    0
29  -0.38612519594749       5   14.00    3    0    0    0    0
30  -0.18817142332233       6   11.00    3    0    0    0    0
31   0.15977647596482       6   14.00    3    0    0    0    0
32   0.53985518513856e-1    7    8.00    3    0    0    0    0
33  -0.28953417958014e-1    7   14.00    3    0    0    0    0
34  -0.13025413381384e-1    8    6.00    3    0    0    0    0
35   0.28948696775778e-2    9    7.00    3    0    0    0    0
36  -0.22647134304796e-2    5   24.00    4    0    0    0    0
37   0.17616456196368e-2    6   22.00    4    0    0    0    0
38   0.58552454482774e-2    2    3.00    0   20  250 1.11    1
39  -0.69251908270028       1    1.00    0   20  375 1.14    1
40   0.15315490030516e+1    2    0.00    0   20  300 1.17    1
41  -0.27380447449783e-2    3    0.00    0   20  225 1.11    1
 ];

FR_N(1:41,Ar)       = Table_30(:,2);
FR_d(1:41,Ar)       = Table_30(:,3);
FR_t(1:41,Ar)       = Table_30(:,4);
FR_c(1:41,Ar)       = Table_30(:,5);
FR_eta(1:41,Ar)     = Table_30(:,6);
FR_beta(1:41,Ar)    = Table_30(:,7);
FR_gamma(1:41,Ar)   = Table_30(:,8);
FR_epsilon(1:41,Ar) = Table_30(:,9);
