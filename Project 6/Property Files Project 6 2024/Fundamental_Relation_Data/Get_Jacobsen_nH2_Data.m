% Loads data into an array for dealing with multiple species.
% Data file for normal hydrogen properties using the data from 
% "Thermodynamic Properties of Cyrogenic Fluids" by R.T. Jacobsen,
% S.G. Penoncello, and E.W. Lemmon, Plenum Press, 1997.
% C.F. Edwards, 1/23/08 

% Be careful:  There is a typo in the book.  The gamma coefficients are not
% listed correctly.  They should only be nonzero where the density exponent
% (l_k) is equal to two.  That is fixed here.  Note that this error is
% embedded in the ALLPROPS data files too--so if you cut & paste data from
% there you will be wrong also.  (This was painful to debug!)

% Constants and ideal-gas ref. state properties:
Ru          = 8314.34;                          % J/kmol-K
M_i(nH2)    = 2.01594;                          % kg/kmol
R_i(nH2)    = Ru/M_i(nH2);                      % J/kg-K
wP_i(nH2)   = -0.214;                           % Pitzer's acentric factor
Tref_i(nH2) = 273.15;                           % K
Pref_i(nH2) = 0.001e6;                          % Pa
rref_i(nH2) = Pref_i(nH2)/R_i(nH2)/Tref_i(nH2); % kg/m3
href = 7747088/M_i(nH2);                        % J/kg
sref = 177838.2/M_i(nH2);                       % J/kg-K

% Fixed-point properties and limits:
Tcrit_i(nH2) = 33.19;               % K
Pcrit_i(nH2) = 1.3152e6;            % Pa
rcrit_i(nH2) = 14.936*M_i(nH2);     % kg/m3

Ttrip_i(nH2)  = 13.95;              % K
Ptrip_i(nH2)  = 0.007199e6;         % Pa
rftrip_i(nH2) = 76.897;             % kg/m3
rgtrip_i(nH2) = 0.13723;            % kg/m3

Tupper_i(nH2) = 500;                % K
Tlower_i(nH2) = 13.95;              % K
Pupper_i(nH2) = 40e6;               % Pa
rupper_i(nH2) = 38.1446*M_i(nH2);   % kg/m3

% Ideal gas Cp coefs. from Table 5.46
Table_5_46 = [
1   -7      0.12155215170e11
2   -6     -0.3639676270e10
3   -5      0.433752654e9
4   -4     -0.2308581738e8
5   -3     -0.386809271e4
6   -2      0.8824013566e5
7   -1     -0.7858708525e4
8    0      0.7248020909e3
9   0.5    -0.1842680629e3
10   1      0.218015504e2
11  1.5    -1.305182
12   2      0.2100317522e-1
13  2.5     0.2391160428e-2
14   3     -0.1824054653e-3
15  3.5     0.5614956073e-5
16   4     -0.7380331013e-7
17   5      0.663577552e-11
    ];

% Translating these into the format of the a0 function gives:
FR_Npoly0(nH2) = 17;
FR_Neinst(nH2) = 0;
FR_Nspec0(nH2) = 0;
FR_t0(1,nH2)    = 0;
FR_t0(2,nH2)    = 1;
FR_t0(3:8,nH2)  = -Table_5_46(1:6,2);
FR_t0(9:17,nH2) = -Table_5_46(9:17,2);
ik = Table_5_46(:,2);
Nk = Table_5_46(:,3);
for(i=1:1:length(Nk))
    ck(i) = -Nk(i)*(Tcrit_i(nH2)^ik(i)) / ( ik(i)*(ik(i)+1) );
end
FR_N0(3:8,nH2)  = ck(1:6);
FR_N0(9:17,nH2) = ck(9:17);
FR_N0(18,nH2)   = Table_5_46(8,3)-1;
FR_N0(19,nH2)   = -Table_5_46(7,3)/Tcrit_i(nH2);

% The constants associated with the zero and first power terms in tau
% set the reference state energy and entropy.  These must be adjusted to
% give the desired state.
% FR_N0(2,nH2) = 0;    % Use this to set h0ref first.
% FR_N0(1,nH2) = 0;    % Then use this to set s0ref (depends on above).
FR_N0(2,nH2) = -444.209519;     % Use this to set h0ref first.
FR_N0(1,nH2) = -633.764269;     % Then use this to set s0ref (depends on above).

% Number of terms of each type in the residual part of the FR:
FR_Npoly(nH2)  = 22; 
FR_Nexp(nH2)   = 0;
FR_Ngaus(nH2)  = 40-22;

% Fundamental Relation coefficients from Table 5.47 of the book:
Table_5_47 = [
1    -8.19267257565         0   3   0   0
2     0.130546404892        0   4   0   0
3    -0.11540517276e-1      0   5   0   0
4     0.839918647516e-1     1   0   0   0
5     1.33747820573         1  0.5  0   0
6    -2.79507010112         1   1   0   0
7     0.482999102128        1   2   0   0
8    -0.148739103766        1   3   0   0
9     0.25598353866e-1      2   0   0   0
10   -0.541546874713e-1     2   1   0   0
11    0.37223474567         2   2   0   0
12    1.8938236094          2   3   0   0
13    0.267165188522e-2     3   0   0   0
14    0.116573604971        3   1   0   0
15   -0.273795055848        3   2   0   0
16   -0.103129506212e-1     4   1   0   0
17    0.397065461542e-1     5   2   0   0
18   -0.840274234789e-1     5   3   0   0
19   -0.110411368495e-1     6   2   0   0
20    0.104800406095e-2     7   2   0   0
21    0.902088706725e-2     7   3   0   0
22   -0.129256450195e-2     8   3   0   0
23    8.19267257565         0   3   2   0.91464479
24   -0.130546404892        0   4   2   0.91464479
25    0.11540517276e-1      0   5   2   0.91464479
26    5.66493981456         2   3   2   0.91464479
27   -0.297031070059        2   4   2   0.91464479
28    0.105554740466e-1     2   5   2   0.91464479
29    1.80209046818         4   3   2   0.91464479
30   -0.135838960943        4   4   2   0.91464479
31    0.394549979521e-1     4   5   2   0.91464479
32    0.353263983498        6   3   2   0.91464479
33   -0.192814451572e-1     6   4   2   0.91464479
34    0.120291028247e-1     6   5   2   0.91464479
35    0.416535430489e-1     8   3   2   0.91464479
36   -0.440891835846e-2     8   4   2   0.91464479
37   -0.980954445376e-3     8   5   2   0.91464479
38    0.53460864973e-2      10  3   2   0.91464479
39    0.275654386288e-3     10  4   2   0.91464479
40   -0.179444975323e-3     10  5   2   0.91464479
 ];

% Note the we are using the Gaussian terms to implement the gamma-weighted
% exponential terms.  As such, we do not need the 5th column above, and the
% gamma coefficient becomes eta below (with the other parts of the Gaussian
% term set to zero).
FR_N(1:40,nH2)       = Table_5_47(:,2);
FR_d(1:40,nH2)       = Table_5_47(:,3);
FR_t(1:40,nH2)       = Table_5_47(:,4);
FR_eta(1:40,nH2)     = Table_5_47(:,6);
FR_beta(1:40,nH2)    = zeros(40,1);
FR_gamma(1:40,nH2)   = zeros(40,1);
FR_epsilon(1:40,nH2) = zeros(40,1);
