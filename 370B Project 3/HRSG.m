function [Tpinch,waterStates,mixStates] = HRSG(mix,waterInput,mdot_m,mdot_w,...
    eff,pr_econ,pr_boil,pr_sup)
% HRSG - Heat Recovery Steam Generator - counterflow heat exchanger
% Assuming no pressure drop on the gas side
% Mining energy from a high temperature gas stream
% Here we don't worry about which is the high/low capacity stream
%
% [Tpinch,waterStates,mixStates] = HRSG(mix,waterInput,mdot_m,mdot_w,eff,pr_econ,pr_boil,pr_sup)
%
% The inputs are as follows
% mix -> GRI30 gas mixture object correctly initialized to gas inlet
% waterInput -> Water liquid-vapor object correctly initialized to water inlet
% mdot_m -> mass flow rate of gas mixture
% mdot_w -> mass flow rate of water
% eff -> overall HRSG effectiveness
% pr_econ -> pressure ratio across economizer (water)
% pr_boil -> pressure ratio across boiler (water)
% pr_sup -> pressure ratio across superheater (water) 
%
% The ouputs are as follows
% Tpinch -> pinch temperature (minimum temp at any point)
% waterStates -> structure containing state parameters of water at all locations
% mixStates -> structure containing state parameters of gas mixture at all locations
%
% w1 -> water inlet
% w2 -> water after economizer (along water stream)
% w3 -> water after boiler (along water stream)
% w4 -> water outlet (along water stream after the superheater)
%
% m1 -> gas mix inlet
% m2 -> gas mix after superheater (along gas stream)
% m3 -> gas mix after boiler (along gas stream)
% m4 -> gas outlet (along gas stream after the economizer)

% get the gas inlet properties
Tm1 = temperature(mix);
Pm1 = pressure(mix);
nm1 = moleFractions(mix);
hm1 = enthalpy_mass(mix);
sm1 = entropy_mass(mix);
xm1 = exergy_mass(mix);
xfm1 = flowExergy_mass(mix);
% initialize a placeholder object for the gas mixture
gasPlaceholder = GRI30;
set(gasPlaceholder,'T',Tm1,'P',Pm1,'X',nm1);

% get the water inlet properties
Tw1 = temperature(waterInput);
Pw1 = pressure(waterInput);
hw1 = enthalpy_mass(waterInput);
sw1 = entropy_mass(waterInput);
xw1 = exergy_mass(waterInput);
xfw1 = flowExergy_mass(waterInput);
% initialize a placeholder object for water
waterPlaceholder = Water;

% get the pressure of water at various stages
Pw2 = Pw1*pr_econ;
Pw3 = Pw2*pr_boil;
Pw4 = Pw3*pr_sup;

% get the pressure of the gas mixture at various stages
% no pressure drop assumption
Pm2 = Pm1;
Pm3 = Pm2;
Pm4 = Pm3;

% find out the maximum possible heat transfer
set(waterPlaceholder,'T',Tm1,'P',Pw4);
hw4_hypothetical = enthalpy_mass(waterPlaceholder);

Qmax = mdot_w*(hw4_hypothetical - hw1);

% find out the actual heat transfer
Qactual = eff*Qmax;

% find out the actual state 4 of water and state 1 of the gas mix
hw4 = hw1 + eff*(hw4_hypothetical - hw1);
set(waterPlaceholder,'H',hw4,'P',Pw4);
Tw4 = temperature(waterPlaceholder);
sw4 = entropy_mass(waterPlaceholder);
xw4 = exergy_mass(waterPlaceholder);
xfw4 = flowExergy_mass(waterPlaceholder);
hm4 = hm1 - Qactual/mdot_m;
set(gasPlaceholder,'H',hm4,'P',Pm4);
equilibrate(gasPlaceholder,'HP');
Tm4 = temperature(gasPlaceholder);
sm4 = entropy_mass(gasPlaceholder);
nm4 = moleFractions(gasPlaceholder);
xm4 = exergy_mass(gasPlaceholder);
xfm4 = flowExergy_mass(gasPlaceholder);

%%% work out the intermediate states
% after water crosses the economizer
set(waterPlaceholder,'P',Pw2,'Vapor',0); % set it to saturated liquid
Tw2 = temperature(waterPlaceholder);
hw2 = enthalpy_mass(waterPlaceholder);
sw2 = entropy_mass(waterPlaceholder);
xw2 = exergy_mass(waterPlaceholder);
xfw2 = flowExergy_mass(waterPlaceholder);
% before water enters the superheater --> saturated liquid
set(waterPlaceholder,'P',Pw3,'Vapor',1);
Tw3 = temperature(waterPlaceholder);
hw3 = enthalpy_mass(waterPlaceholder);
sw3 = entropy_mass(waterPlaceholder);
xw3 = exergy_mass(waterPlaceholder);
xfw3 = flowExergy_mass(waterPlaceholder);
% after the gas crosses the superheater
Q_sup = mdot_w*(hw4 - hw3);
hm2 = hm1 - Q_sup/mdot_m;
set(gasPlaceholder,'H',hm2,'P',Pm2);
equilibrate(gasPlaceholder,'HP');
Tm2 = temperature(gasPlaceholder);
sm2 = entropy_mass(gasPlaceholder);
nm2 = moleFractions(gasPlaceholder);
xm2 = exergy_mass(gasPlaceholder);
xfm2 = flowExergy_mass(gasPlaceholder);
% after the gas crosses the boiler
Q_boil = mdot_w*(hw3 - hw2);
hm3 = hm2 - Q_boil/mdot_m;
set(gasPlaceholder,'H',hm3,'P',Pm2);
equilibrate(gasPlaceholder,'HP');
Tm3 = temperature(gasPlaceholder);
sm3 = entropy_mass(gasPlaceholder);
nm3 = moleFractions(gasPlaceholder);
xm3 = exergy_mass(gasPlaceholder);
xfm3 = flowExergy_mass(gasPlaceholder);

% calculate temperature differences at all locations
% Water outlet, after water boiler, after water economizer, water inlet
TpinchAllLocs = [Tm1-Tw4;Tm2-Tw3;Tm3-Tw2;Tm4-Tw1];


% return value of the pinch temperature
Tpinch = min(TpinchAllLocs);


% return the states of water at various locations
waterStates.w1.location = 'Water inlet';
waterStates.w1.T = Tw1;
waterStates.w1.P = Pw1;
waterStates.w1.h = hw1;
waterStates.w1.s = sw1;
waterStates.w1.x = xw1;
waterStates.w1.xf = xfw1;

waterStates.w2.location = 'After economizer';
waterStates.w2.T = Tw2;
waterStates.w2.P = Pw2;
waterStates.w2.h = hw2;
waterStates.w2.s = sw2;
waterStates.w2.x = xw2;
waterStates.w2.xf = xfw2;

waterStates.w3.location = 'After boiler';
waterStates.w3.T = Tw3;
waterStates.w3.P = Pw3;
waterStates.w3.h = hw3;
waterStates.w3.s = sw3;
waterStates.w3.x = xw3;
waterStates.w3.xf = xfw3;

waterStates.w4.location = 'Water outlet';
waterStates.w4.T = Tw4;
waterStates.w4.P = Pw4;
waterStates.w4.h = hw4;
waterStates.w4.s = sw4;
waterStates.w4.x = xw4;
waterStates.w4.xf = xfw4;

% return the states of the gas at various locations
mixStates.m1.location = 'Gas mix inlet';
mixStates.m1.T = Tm1;
mixStates.m1.P = Pm1;
mixStates.m1.h = hm1;
mixStates.m1.s = sm1;
mixStates.m1.n = nm1;
mixStates.m1.x = xm1;
mixStates.m1.xf = xfm1;

mixStates.m2.location = 'Gas mix after superheater';
mixStates.m2.T = Tm2;
mixStates.m2.P = Pm2;
mixStates.m2.h = hm2;
mixStates.m2.s = sm2;
mixStates.m2.n = nm2;
mixStates.m2.x = xm2;
mixStates.m2.xf = xfm2;

mixStates.m3.location = 'Gas mix after boiler';
mixStates.m3.T = Tm3;
mixStates.m3.P = Pm3;
mixStates.m3.h = hm3;
mixStates.m3.s = sm3;
mixStates.m3.n = nm3;
mixStates.m3.x = xm3;
mixStates.m3.xf = xfm3;

mixStates.m4.location = 'Gas mix outlet';
mixStates.m4.T = Tm4;
mixStates.m4.P = Pm4;
mixStates.m4.h = hm4;
mixStates.m4.s = sm4;
mixStates.m4.n = nm4;
mixStates.m4.x = xm4;
mixStates.m4.xf = xfm4;

end