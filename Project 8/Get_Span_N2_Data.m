% Loads data into an array for dealing with multiple species.
% Data file for nitrogen properties using the data from 
% "Thermodynamic Properties of Nitrogen" by Span et al.,
% J. Phy. Chem. Ref. Data, 29(6), 1361, 2000
% C.F. Edwards, 9/2/07 

% Constants and ideal-gas ref. state properties:
Ru         = 8314.51;                       % J/kmol-K
M_i(N2)    = 28.01348;                      % kg/kmol
R_i(N2)    = Ru/M_i(N2);                    % J/kg-K
wP_i(N2)   = 0.037;                         % Pitzer's acentric factor
Tref_i(N2) = 298.15;                        % K
Pref_i(N2) = 101325;                        % Pa
rref_i(N2) = Pref_i(N2)/R_i(N2)/Tref_i(N2); % kg/m3

% Fixed-point properties and limits:
Tcrit_i(N2) = 126.192;                      % K
Pcrit_i(N2) = 3.3958e6;                     % Pa
rcrit_i(N2) = 11.1839*M_i(N2);              % kg/m3

Ttrip_i(N2)  = 63.151;                      % K
Ptrip_i(N2)  = 12523;                       % Pa
rftrip_i(N2) = 30.9573*M_i(N2);             % kg/m3
rgtrip_i(N2) = 0.02407*M_i(N2);             % kg/m3

Tupper_i(N2) = 1000;                        % K
Tlower_i(N2) = 63.151;                      % K
Pupper_i(N2) = 2200e6;                      % Pa
rupper_i(N2) = 41.528*M_i(N2);              % kg/m3

% Ideal gas Helmholtz fit coefficients for Eq.4.6.
FR_Npoly0(N2) = 5;
FR_Neinst(N2) = 1;
FR_t0(1:5,N2) = [0; 1; -1; -2; -3];
% The original Span et al. paper uses nonzero values for h0ref and s0ref.
% Use these and the associated FR_N0 values if you want to compare with the
% tables in the original paper.  The Lemmon air model uses these same
% values.
href_i(N2) = 8670000./M_i(N2);              % J/kg
sref_i(N2) = 191500./M_i(N2);               % J/kg-K
% Note that the order of the coefficients is a little different from that
% in the Span paper:  We put the log term after the polynomial terms and
% before the Einstein term(s).
FR_N0(1:8,N2) = [
   -12.76952708
    -0.00784163
    -1.934819e-4
    -1.247742e-5
     6.678326e-8
    2.5
    0.0
    1.012941
    ];
% The gamma coefficient is listed as a8 in the original paper.
FR_gamma0(8,N2) = 26.65788;

% Number of terms of each type in the residual part of the FR:
FR_Npoly(N2)  = 6; 
FR_Nexp(N2)   = 32-6;
FR_Ngaus(N2)  = 36-32;

% Fundamental Relation coefficients from Table 17 of the paper:
Table_17 = [
1    0.924803575275       1    0.25    0    0    0    0    0
2   -0.492448489428       1    0.875   0    0    0    0    0
3    0.661883336938       2    0.5     0    0    0    0    0
4   -0.192902649210e+1    2    0.875   0    0    0    0    0
5   -0.622469309629e-1    3    0.375   0    0    0    0    0
6    0.349943957581       3    0.75    0    0    0    0    0
7    0.564857472498       1    0.5     1    0    0    0    0
8   -0.161720005987e+1    1    0.75    1    0    0    0    0
9   -0.481395031883       1    2.00    1    0    0    0    0
10   0.421150636384       3    1.25    1    0    0    0    0
11  -0.161962230825e-1    3    3.5     1    0    0    0    0
12   0.172100994165       4    1.0     1    0    0    0    0
13   0.735448924933e-2    6    0.5     1    0    0    0    0
14   0.168077305479e-1    6    3.0     1    0    0    0    0
15  -0.107626664179e-2    7    0.0     1    0    0    0    0
16  -0.137318088513e-1    7    2.75    1    0    0    0    0
17   0.635466899859e-3    8    0.75    1    0    0    0    0
18   0.304432279419e-2    8    2.5     1    0    0    0    0
19  -0.435762336045e-1    1    4.0     2    0    0    0    0
20  -0.723174889316e-1    2    6.0     2    0    0    0    0
21   0.389644315272e-1    3    6.0     2    0    0    0    0
22  -0.212201363910e-1    4    3.0     2    0    0    0    0
23   0.408822981509e-2    5    3.0     2    0    0    0    0
24  -0.551990017984e-4    8    6.0     2    0    0    0    0
25  -0.462016716479e-1    4   16.0     3    0    0    0    0
26  -0.300311716011e-2    5   11.0     3    0    0    0    0
27   0.368825891208e-1    5   15.0     3    0    0    0    0
28  -0.255856846220e-2    8   12.0     3    0    0    0    0
29   0.896915264558e-2    3   12.0     4    0    0    0    0
30  -0.441513370350e-2    5    7.0     4    0    0    0    0
31   0.133722924858e-2    6    4.0     4    0    0    0    0
32   0.264832491957e-3    9   16.0     4    0    0    0    0
33   0.196688194015e+2    1    0.0     2   20  325 1.16    1
34  -0.209115600730e+2    1    1.0     2   20  325 1.16    1
35   0.167788306989e-1    3    2.0     2   15  300 1.13    1
36   0.262767566274e+4    2    3.0     2   25  275 1.25    1
 ];

FR_N(1:36,N2)       = Table_17(:,2);
FR_d(1:36,N2)       = Table_17(:,3);
FR_t(1:36,N2)       = Table_17(:,4);
FR_c(1:36,N2)       = Table_17(:,5);
FR_eta(1:36,N2)     = Table_17(:,6);
FR_beta(1:36,N2)    = Table_17(:,7);
FR_gamma(1:36,N2)   = Table_17(:,8);
FR_epsilon(1:36,N2) = Table_17(:,9);
