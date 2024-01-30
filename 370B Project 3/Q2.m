%% Close Everything
clc;clear;close all;
fprintf('\nClosing everything...\n')

%% Declare Global Variables
global To Po xo mu_o g_tm;                      % Formatted according to Chris' LHV/HHV scripts
g_tm=0;
To=293.15;
Po=101325;

%% Declare Local Variables and Assign Values
fprintf('\nDeclaring Local Variables...\n')
TIT=1700;                                       % Turbine Inlet Temperature from Q2
PR=23;                                          % Pressure Ratio from Q2
BurnerPR= 0.95;                                 % Values from Project 1
AirComp_eff=0.90;
FuelComp_eff=0.90;
Turbine_eff=0.88;
Air=zeros(1, 4);                                % State variables order-> Enthalpy, Entropy, specific exergy, flow exergy
Fuel=zeros(1, 4);                               % Enthalpy, Entropy, Exergy (Not Specific), Flow Exergy (Not Specific)
Mix=zeros(1, 4);                                % Same as above

%% Declare Cantera Object
gas=Solution('gri30.yaml');
N=nSpecies(gas);
% Gases-> N2, O2, CO2, AR, H20, CH4, C2H6, C3H8
indx=[speciesIndex(gas,'N2'), speciesIndex(gas,'O2'), speciesIndex(gas,'CO2'), speciesIndex(gas,'AR'), speciesIndex(gas,'H2O'), speciesIndex(gas,'CH4'), speciesIndex(gas,'C2H6'), speciesIndex(gas,'C3H8')];
%% Declare Fuel (Copied from Chris's Project 1 script)
xfuel=zeros(1,N);
xfuel(indx(1:8))=[0.018, 0.010, 0.010, 0, 0, 0.907, 0.036, 0.019];
%% Calculate Fuel LHV
fprintf('\nCalculating LHV for specified Fuel...\n')
set(gas,'T',To,'P',Po,'X',xfuel);
Fuel_LHV=LHV_mass(gas);
%% Declare Dry Air 
xo=zeros(1,N);
xo(indx(1:4))=[0.78084, 0.2094, 0.000412, 0.0093];% Dry Air composition
for i=1:1:N                                     % Re-normalize composition
    xo(i)=xo(i)/sum(xo);
end
%% Add humidity to the air
w=Water;
setState_Tsat(w,[To 1]);
N_new=1/(1-0.5*pressure(w)/Po);                  % New total N after RH=50% water is added
for i=1:1:N                                      % Re-normalize to add humidity
    xo(i)=xo(i)/N_new;
end
xo(indx(5))=0.5*pressure(w)/Po;                  % Add humidity
set(gas,'T',To,'P',Po,'X',xo);
mu_o=chemPotentials(gas);
%% Stations (According to Project 1 Schematic)

% Station Air-1
fprintf('\nStation Air-1...\n')
set(gas,'T',To,'P',Po,'X',xo);
MMW_air=meanMolecularWeight(gas);
mdot_air=1; % Set to 1kg/s
Air(1, 1:4)=[enthalpy_mass(gas), entropy_mass(gas), exergy_mass(gas)*mdot_air, flowExergy_mass(gas)*mdot_air];

%Station Air-2 Polytropic Compression (Copied from Chris' Project 1 Script)
fprintf('\nCompressing Air...\n')
P2=Po*PR;
dP=(P2-Po)/50;
s=entropy_mass(gas);
h=enthalpy_mass(gas);
for P=Po:dP:P2
    set(gas,'S',s,'P',P);
    hs=enthalpy_mass(gas);
    dhs=hs-h; 
    h=dhs/AirComp_eff+h;
    set(gas,'H',h,'P',P);
    s=entropy_mass(gas);
end
Air(2, 1:4)=[enthalpy_mass(gas), entropy_mass(gas), exergy_mass(gas)*mdot_air, flowExergy_mass(gas)*mdot_air];

% Station Fuel-1
fprintf('\nStation Fuel-1...\n')
set(gas,'T',To,'P',Po,'X',xfuel);                   % Start at ambient
MMW_fuel=meanMolecularWeight(gas);
Fuel(1, 1:4)=[enthalpy_mass(gas), entropy_mass(gas), exergy_mass(gas), flowExergy_mass(gas)];

% Station Fuel-2
% Polytropic Compression (Copied from Chris' Project 1 Script)
fprintf('\nCompressing Fuel...\n')
Pf2=Po*2*PR;
dP=(Pf2-Po)/50;
s=Fuel(1, 2);
h=Fuel(1, 1);
for P=Po:dP:Pf2
    set(gas,'S',s,'P',P);
    hs=enthalpy_mass(gas);
    dhs=hs-h; 
    h=dhs/FuelComp_eff+h;
    set(gas,'H',h,'P',P);
    s=entropy_mass(gas);
end
Fuel(2, 1:4)=[enthalpy_mass(gas), entropy_mass(gas), exergy_mass(gas), flowExergy_mass(gas)];

% TIT Validation for correct mdot_fuel (Copied from Chris' Project 1 Script)
fprintf('\nAdjusting Fuel Flow Rate till correct inlet temperature is found...\n')
dfuel=(mdot_air/15)/2000;
for mdot_fuel=dfuel:dfuel:mdot_air
    mdot_mix=mdot_air+mdot_fuel;
    h=(mdot_fuel*Fuel(2, 1)+mdot_air*Air(2, 1))/mdot_mix;
    xmix=(mdot_fuel/MMW_fuel)*xfuel+(mdot_air/MMW_air)*xo;
    set(gas,'H',h,'P',P2*BurnerPR,'X',xmix);
    Mix(1, 3)=exergy_mass(gas)*mdot_mix;
    Mix(1, 4)=flowExergy_mass(gas)*mdot_mix;
    set(gas,'H',h,'P',P2*BurnerPR);
    equilibrate(gas,'HP');
    if(temperature(gas)>=TIT)
        fprintf('\nTurbine Inlet Temperature Matched!...\n')
        break
    end
end
fprintf('\nStation Mix-2...\n')
Mix(2, 1:4)=[enthalpy_mass(gas), entropy_mass(gas), exergy_mass(gas)*mdot_mix, flowExergy_mass(gas)*mdot_mix];

% Convert Fuel specific exergies to total exergies
Fuel(1:2,3)=Fuel(1:2, 3).*mdot_fuel;
Fuel(1:2,4)=Fuel(1:2, 4).*mdot_fuel;

% Station Mix-3
% Polytropic Expansion (Copied from Chris' Project 1 Script)
fprintf('\nStation Mix-3...\n')
dP=(P2*BurnerPR-Po)/50;
h=Mix(2, 1);
s=Mix(2, 2);
for P=P2*BurnerPR:-dP:Po
    set(gas,'S',s,'P',P);
    hs=enthalpy_mass(gas);
    dhs=h-hs; 
    h=h-dhs*Turbine_eff;
    set(gas,'H',h,'P',P);
    s=entropy_mass(gas);
end
Mix(3, 1:4)=[enthalpy_mass(gas), entropy_mass(gas), exergy_mass(gas)*mdot_mix, flowExergy_mass(gas)*mdot_mix];

%% Work Calculations
% Work-> Air compressor, Fuel Compressor, Turbine, Overall, 
work=[mdot_air*(Air(2, 1)-Air(1, 1)), mdot_fuel*(Fuel(2, 1)-Fuel(1, 1)), mdot_mix*(Mix(2, 1)-Mix(3, 1))];         % Air Compressor
work(4)=work(3)-work(2)-work(1);                % Overall Work
Eff_LHV=work(4)/(mdot_fuel*Fuel_LHV);           % First-Law Efficiency
mdot_fuel=mdot_fuel*1e6;
%% Plotting
% Exergy Path
fprintf('\nPlotting Exergy Paths...\n')
figure(1)
hold on
plot([1,2,3],[Air(1:2, 3)', Mix(1, 3)]/mdot_fuel,'-ob');
plot([1,2,3],[Fuel(1:2, 3)', Mix(1, 3)]/mdot_fuel,'-or');
plot([3,4,5],Mix(1:3, 3)/mdot_fuel,'-ok');
xlabel('Station');
ylabel('Fuel-Specific Exergy MJ/Kg-Fuel');
legend('Air','Fuel','Mix');
plotfixer
hold off

% Flow Exergy Path
figure(2)
hold on
plot([1,2,3],[Air(1:2, 4)', Mix(1, 4)]/mdot_fuel,'-ob');
plot([1,2,3],[Fuel(1:2, 4)', Mix(1, 4)]/mdot_fuel,'-or');
plot([3,4,5],Mix(1:3, 4)/mdot_fuel,'-ok');
xlabel('Station');
ylabel('Fuel-Specific Flow Exergy MJ/Kg-Fuel');
legend('Air','Fuel','Mix');
hold off
plotfixer

%% Exergy Balance Plot
fprintf('\nPerforming Exergy Balance...\n')
data=[work(2)+Fuel(1,4)-Fuel(2,4), work(1)+Air(1, 4)-Air(2, 4), Fuel(2,4)+Air(2, 4)-Mix(1, 4), Mix(1, 4)-Mix(2, 4), Mix(2, 4)-work(3)-Mix(3, 4), Mix(3, 4)];
figure(3)
bar(100*data/sum(data));
xlabel('Conversion Device');
ylabel('Eexergy Loss(%)');
text(.1,45,sprintf('Efficiency (LHV):%f',100*Eff_LHV));
text(.1,40,sprintf('Exergy Efficiency:%f',100*work(4)/(Fuel(1,4)+Air(1, 4))));
text(0.1,35,'Device Key:');
text(0.1,30,'1-Fuel Compressor');
text(0.1,25,'2-Air Compressor');
text(0.1,20,'3-Premixer');
text(0.1,15,'4-Combustor');
text(0.1,10,'5-Turbine');
text(0.1,5,'6-Exhaust');
plotfixer
legend('off');