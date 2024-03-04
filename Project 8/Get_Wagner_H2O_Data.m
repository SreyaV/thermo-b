% Loads data into an array for dealing with multiple species.
% Data file for IAPWS-95 water properties using the data from 
% Wagner and Pruss, "Thermodynamic Properties of
% Ordinary Water," J. Phys. Chem. Ref. Data, 31(2), 387 ,2002
% Adam Berger, 7/10/2007
% Edited by C.F. Edwards, 9/9/2007

% Constants and ideal-gas ref. state properties:
M_i(H2O)    = 18.015268;                  % kg/kmol
Ru          = 8314.472;                   % J/kmol-K
R_i(H2O)    = 461.51805;                  % J/kg-K
wP_i(H2O)   = 0.3443;                     % Pitzer's acentric factor
% Wagner and Pruss do not use an ideal gas reference state.  See below.

% Fixed-point properties and limits:
Tcrit_i(H2O) = 647.096;                   % K
Pcrit_i(H2O) = 22.064e6;                  % Pa
rcrit_i(H2O) = 322;                       % kg/m3

Ttrip_i(H2O)  = 273.16;                   % K
Ptrip_i(H2O)  = 611.655;                  % Pa
rftrip_i(H2O) = 999.793;                  % kg/m3
rgtrip_i(H2O) = .00485458;                % kg/m3

Tupper_i(H2O) = 1273;                     % K
Tlower_i(H2O) = 251.165;                  % K
Pupper_i(H2O) = 1000e6;                   % Pa
rupper_i(H2O) = 1237.39;                  % kg/m3

% Ideal gas Helmholtz fit coefficients for Eq.4.6.
FR_Npoly0(H2O) = 2;
FR_Neinst(H2O) = 5;
FR_t0(1:2,H2O) = [0; 1];
% The original Wagner and Pruss paper uses zero values for href and sref at
% the saturated liquid state on the triple line.  Note that this is NOT the
% ideal gas state, but the real fluid state that is set to zero.
FR_N0(1:9,H2O) = [
    -8.32044648201
     6.6832105268
     3.00632
     0.0
     0.012436
     0.97315
     1.27950
     0.96956
     0.24873
    ];
% The first two coefficients of FR_N0 determine the reference state for
% enthalpy and entropy.  The values given above will give
% values that are nominally zero for u and s at the sat liquid state on
% the triple line.  To get these to be exactly zero, you must 
% adjust the constants to get zero to the precision required.
% Adjust coefficient 2 to get uf correct at TP:
FR_N0(2,H2O) = FR_N0(2,H2O) + 7.6385e-10;   % Gets uf = 0
% Adjust coefficient 1 to get sf correct at TP:
FR_N0(1,H2O) = FR_N0(1,H2O) - 1.67e-9;      % Gets sf = 0
% The gamma coefficients are below:
FR_gamma0(5:9,H2O) = [
    1.28728967
    3.53734222
    7.74073708
    9.24437796
   27.5075105
    ];

% Number of terms of each type in the residual part of the FR:
FR_Npoly(H2O)  = 7; 
FR_Nexp(H2O)   = 51-7;
FR_Ngaus(H2O)  = 54-51;
FR_Nnonan(H2O) = 56-54;

% Fundamental Relation coefficients entered from Wagner and Pruss Table 6.2:
DATA = [
1    0.12533547935523e-1    0    1    -0.5      %Polynomial terms
2    0.78957634722828e1     0    1    0.875  
3   -0.87803203303561e1     0    1    1         
4    0.31802509345418       0    2    0.5    
5   -0.26145533859358       0    2    0.75   
6   -0.78199751687981e-2    0    3    0.375  
7    0.88089493102134e-2    0    4    1      
 
8   -0.66856572307965       1    1    4         %Exponential terms
9    0.20433810950965       1    1    6      
10  -0.66212605039687e-4    1    1    12     
11  -0.19232721156002       1    2    1      
12  -0.25709043003438       1    2    5      
13   0.16074868486251       1    3    4      
14  -0.40092828925807e-1    1    4    2      
15   0.39343422603254e-6    1    4    13     
16  -0.75941377088144e-5    1    5    9      
17   0.56250979351888e-3    1    7    3      
18  -0.15608652257135e-4    1    9    4      
19   0.11537996422951e-8    1    10   11     
20   0.36582165144204e-6    1    11   4      
21  -0.13251180074668e-11   1    13   13     
22  -0.62639586912454e-9    1    15   1      
23  -0.10793600908932       2    1    7      
24   0.17611491008752e-1    2    2    1      
25   0.22132295167546       2    2    9      
26  -0.40247669763528       2    2    10     
27   0.58083399985759       2    3    10     
28   0.49969146990806e-2    2    4    3      
29  -0.31358700712549e-1    2    4    7      
30  -0.74315929710341       2    4    10     
31   0.47807329915480       2    5    10     
32   0.20527940895948e-1    2    6    6      
33  -0.13636435110343       2    6    10     
34   0.14180634400617e-1    2    7    10     
35   0.83326504880713e-2    2    9    1      
36  -0.29052336009585e-1    2    9    2      
37   0.38615085574206e-1    2    9    3      
38  -0.20393486513704e-1    2    9    4      
39  -0.16554050063734e-2    2    9    8      
40   0.19955571979541e-2    2    10   6      
41   0.15870308324157e-3    2    10   9      
42  -0.16388568342530e-4    2    12   8      
43   0.43613615723811e-1    3    3    16    
44   0.34994005463765e-1    3    4    22     
45  -0.76788197844621e-1    3    4    23     
46   0.22446277332006e-1    3    5    23     
47  -0.62689710414685e-4    4    14   10     
48  -0.55711118565645e-9    6    3    50     
49  -0.19905718354408       6    6    44     
50   0.31777497330738       6    6    46     
51  -0.11841182425981       6    6    50      

52  -0.31306260323435e2     0    3    0         %Gaussian terms 
53   0.31546140237781e2     0    3    1      
54  -0.25213154341695e4     0    3    4      

55  -0.14874640856724       0    0    0         %Nonanalytic Terms
56   0.31806110878444       0    0    0      
];

DATA_gaus = [
52  20  150   1.21    1
53  20  150   1.21    1
54  20  250   1.25    1
];

DATA_nonan = [
55  3.5   0.85  .2  28  700   .32   .3  
56  3.5   0.95  .2  32  800   .32   .3
];

FR_N(1:56,H2O)        = DATA(:,2);
FR_d(1:56,H2O)        = DATA(:,4);
FR_t(1:56,H2O)        = DATA(:,5);
FR_c(1:56,H2O)        = DATA(:,3);
FR_eta(52:54,H2O)     = DATA_gaus(:,2);     % This is the alpha coefficient in the original paper.
FR_beta(52:54,H2O)    = DATA_gaus(:,3);
FR_gamma(52:54,H2O)   = DATA_gaus(:,4);
FR_epsilon(52:54,H2O) = DATA_gaus(:,5);

FR_a(55:56,H2O)       = DATA_nonan(:,2);
FR_b(55:56,H2O)       = DATA_nonan(:,3);
FR_NAbeta(55:56,H2O)  = DATA_nonan(:,8);    % NAbeta is the beta coefficient in the Non-Analytic terms.
FR_A(55:56,H2O)       = DATA_nonan(:,7);
FR_B(55:56,H2O)       = DATA_nonan(:,4);
FR_C(55:56,H2O)       = DATA_nonan(:,5);
FR_D(55:56,H2O)       = DATA_nonan(:,6);
