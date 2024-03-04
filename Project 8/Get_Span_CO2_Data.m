% Loads data into an array for dealing with multiple species.
% Data file for carbon dioxide properties using the data from 
% "A New Equation of State for Carbon Dioxide..." by Span, and Wagner,
% J. Phy. Chem. Ref. Data, 25(6), 1509, 1996
% Adam Berger, 6/29/2007 
% Modified by C.F. Edwards, 9/3/2007

% Constants and ideal-gas ref. state properties:
M_i(CO2)    = 44.0098;                      % kg/kmol
Ru          = 8314.51;                      % J/kmol-K
R_i(CO2)    = 188.9241;                     % J/kg-K
wP_i(CO2)   = 0.224;                        % Pitzer's acentric factor 
Tref_i(CO2) = 298.15;                       % K
Pref_i(CO2) = 101325;                       % Pa
rref_i(CO2) = Pref_i(CO2)/R_i(CO2)/Tref_i(CO2); % kg/m3

% Fixed-point properties and limits:
Tcrit_i(CO2) = 304.1282;                    % K
Pcrit_i(CO2) = 7.3773e6;                    % Pa
rcrit_i(CO2) = 467.6;                       % kg/m3

Ttrip_i(CO2)  = 216.592;                    % K
Ptrip_i(CO2)  = 517950;                     % Pa
rftrip_i(CO2) = 1178.46;                    % kg/m3
rgtrip_i(CO2) = 13.761;                     % kg/m3

Tupper_i(CO2) = 1100;                       % K
Tlower_i(CO2) = 186.436;                    % K
Pupper_i(CO2) = 800e6;                      % Pa
rupper_i(CO2) = 1495.70;                    % kg/m3

% Ideal gas Helmholtz fit coefficients for Eq.4.6.
FR_Npoly0(CO2) = 2;
FR_Neinst(CO2) = 5;
FR_t0(1:2,CO2) = [0; 1];
% The original Span and Wagner paper uses zero values for h0ref and s0ref.
% Use these and the associated FR_N0 values if you want to compare with the
% tables in the original paper.
href_i(CO2) = 0;                             % J/kg
sref_i(CO2) = 0;                             % J/kg-K
% Note that the order of the coefficients is a little different from that
% in the paper:  We put the log term after the polynomial terms and
% before the Einstein terms.
FR_N0(1:9,CO2) = [
     8.37304456
    -3.70454304
     2.5
     0.0
     1.99427042
     0.62105248
     0.41195293
     1.04028922
     0.08327678
    ];
% The first two coefficients of FR_N0 determine the reference state for
% enthalpy and entropy.  The values given above will give
% values that are nominally zero for h and s for an ideal gas at T = 298.15 K
% and P = 101325 Pa.  To get these to be exactly zero, you must 
% adjust the constants to get zero to the precision required.
% Adjust coefficient 2 to get h correct:
FR_N0(2,CO2) = FR_N0(2,CO2) + 1.4105e-9;  % Gets h0 = 0 at 1 atm, 25C
% Adjust coefficient 1 to get s correct:
FR_N0(1,CO2) = FR_N0(1,CO2) + 4.85e-9;    % Gets s0 = 0 at 1 atm, 25C
% The gamma coefficients are listed as theta0i in the original paper.
FR_gamma0(5:9,CO2) = [
    3.15163
    6.11190
    6.77708
   11.32384
   27.08792
    ];

% Number of terms of each type in the residual part of the FR:
FR_Npoly(CO2)  = 7; 
FR_Nexp(CO2)   = 34-7;
FR_Ngaus(CO2)  = 39-34;
FR_Nnonan(CO2) = 42-39;

% Fundamental Relation coefficients entered from Span and Wagner Table 31:
Table_31 = [
1    0.38856823203161       1    0.00    0     0     0     0     0
2    0.29385475942740e+1    1    0.75    0     0     0     0     0
3   -0.55867188534934e+1    1    1.00    0     0     0     0     0 
4   -0.76753199592477       1    2.00    0     0     0     0     0
5    0.31729005580416       2    0.75    0     0     0     0     0
6    0.54803315897767       2    2.00    0     0     0     0     0
7    0.12279411220335       3    0.75    0     0     0     0     0
  
8    0.21658961543220e+1    1    1.50    1     0     0     0     0
9    0.15841735109724e+1    2    1.50    1     0     0     0     0
10  -0.23132705405503       4    2.50    1     0     0     0     0
11   0.58116916431436e-1    5    0.00    1     0     0     0     0
12  -0.55369137205382       5    1.50    1     0     0     0     0
13   0.48946615909422       5    2.00    1     0     0     0     0
14  -0.24275739843501e-1    6    0.00    1     0     0     0     0
15   0.62494790501678e-1    6    1.00    1     0     0     0     0
16  -0.12175860225246       6    2.00    1     0     0     0     0
17  -0.37055685270086       1    3.00    2     0     0     0     0
18  -0.16775879700426e-1    1    6.00    2     0     0     0     0
19  -0.11960736637987       4    3.00    2     0     0     0     0
20  -0.45619362508778e-1    4    6.00    2     0     0     0     0
21   0.35612789270346e-1    4    8.00    2     0     0     0     0
22  -0.74427727132052e-2    7    6.00    2     0     0     0     0
23  -0.17395704902432e-2    8    0.00    2     0     0     0     0
24  -0.21810121289527e-1    2    7.00    3     0     0     0     0
25   0.24332166559236e-1    3   12.00    3     0     0     0     0
26  -0.37440133423463e-1    3   16.00    3     0     0     0     0
27   0.14338715756878       5   22.00    4     0     0     0     0
28  -0.13491969083286       5   24.00    4     0     0     0     0
29  -0.23151225053480e-1    6   16.00    4     0     0     0     0
30   0.12363125492901e-1    7   24.00    4     0     0     0     0
31   0.21058321972940e-2    8    8.00    4     0     0     0     0
32  -0.33958519026368e-3   10    2.00    4     0     0     0     0
33   0.55993651771592e-2    4   28.00    5     0     0     0     0
34  -0.30335118055646e-3    8   14.00    6     0     0     0     0
 
35  -0.21365488688320e+3    2    1.00   25   325  1.16  1.00     0
36   0.26641569149272e+5    2    0.00   25   300  1.19  1.00     0
37  -0.24027212204557e+5    2    1.00   25   300  1.19  1.00     0
38  -0.28341603423999e+3    3    3.00   15   275  1.25  1.00     0
39   0.21247284400179e+3    3    3.00   20   275  1.22  1.00     0
  
40  -0.66642276540751     3.5   0.875  0.3   0.7   0.3  10.0   275
41   0.72608632349897     3.5   0.925  0.3   0.7   0.3  10.0   275 
42   0.55068668612842e-1  3.0   0.875  0.3   0.7   1.0  12.5   275
 ];

FR_N(1:42,CO2)        = Table_31(:,2);
FR_d(1:39,CO2)        = Table_31(1:39,3);
FR_t(1:39,CO2)        = Table_31(1:39,4);
FR_c(8:34,CO2)        = Table_31(8:34,5);
FR_eta(35:39,CO2)     = Table_31(35:39,5);  % This is the alpha coefficient in the original paper.
FR_beta(35:39,CO2)    = Table_31(35:39,6);
FR_gamma(35:39,CO2)   = Table_31(35:39,7);
FR_epsilon(35:39,CO2) = Table_31(35:39,8);
FR_a(40:42,CO2)       = Table_31(40:42,3);
FR_b(40:42,CO2)       = Table_31(40:42,4);
FR_NAbeta(40:42,CO2)  = Table_31(40:42,5);  % NAbeta is the beta coefficient in the Non-Analytic terms.
FR_A(40:42,CO2)       = Table_31(40:42,6);
FR_B(40:42,CO2)       = Table_31(40:42,7);
FR_C(40:42,CO2)       = Table_31(40:42,8);
FR_D(40:42,CO2)       = Table_31(40:42,9);
