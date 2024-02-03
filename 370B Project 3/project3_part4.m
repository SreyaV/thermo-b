%% Project 3 Part 4
close all;clear;clc;
format long e;
%% Load gas and steam turbine inputs
TIT = 1700; % Turbine Inlet Temperature (in K)
GPR = 23; % Gas turbine cycle Pressure Ratio --> air compressor side
burner_pratio = 0.95; % burner pressure ratio
eta_airComp = 0.90; % air compressor polytropic efficiency
eta_fuelComp = 0.90; % fuel compressor polytropic efficiency
eta_gasTurb = 0.88; % gas turbine polytropic efficiency
eta_feedPump = 0.85; % water feed pump polytropic efficiency
HRSG_eff = 0.85; % HRSG effectiveness
econ_pratio = 0.92; % economizer pressure ratio
boil_pratio = 0.92; % boiler pressure ratio
superheat_pratio = 0.92; % superheater pressure ratio
mdot_air = 1; % air mass flow rate to gas turbine engine (in kg/s)

%% set environmental parameters
Patm = 1*oneatm; % atmospheric pressure
Tatm = 25+273.15; % atmospheric temperature
RHatm = 0.5; % relative humidity
%% set number of steps for compressors and turbines
nsteps = 250;

%% set air and fuel composition
air = GRI30; % object for air
fuel = GRI30; % object for fuel
gasmix = GRI30; % object for air-fuel mixture
airgasnames = {'N2','O2','AR','CO2','H2O'};
airgasComposition = [0.78084 0.2094 0.0093 0.000412 0];
fuelgasnames = {'CH4','C2H6','C3H8','N2','CO2','O2'};
fuelgasComposition = [0.907 0.036 0.019 0.018 0.01 0.01];

nsp = nSpecies(air);
xair = zeros(nsp,1);
xfuel = xair;

% normalize it
airgasComposition = airgasComposition/sum(airgasComposition);
% add humidity
waterObj = Water;
Psat = satPressure(waterObj,Tatm);
xWater = RHatm*Psat/Patm;
Nnew = 1/(1-xWater);
airgasComposition = airgasComposition/Nnew;
airgasComposition(end) = xWater;
% renormalize
airgasComposition = airgasComposition/sum(airgasComposition);

% set air and fuel compositions
iair = speciesIndex(air,airgasnames);
xair(iair) = airgasComposition;
ifuel = speciesIndex(fuel,fuelgasnames);
xfuel(ifuel) = fuelgasComposition;

set(air,'T',Tatm,'P',Patm,'X',xair);
set(fuel,'T',Tatm,'P',Patm,'X',xfuel);

%% Compress air 
P1a = Patm;
T1a = Tatm;
h1a = enthalpy_mass(air);
P2a = GPR*Patm;

fprintf('Compressing air\n');
[PairComp,TairComp,state_2a,path_1a_2a] = compTurb(air,P1a,T1a,P2a,eta_airComp,nsteps);
T2a = state_2a.T;
h2a = state_2a.H;

%% Compress fuel
P1f = Patm;
T1f = Tatm;
h1f = enthalpy_mass(fuel);
P2f = 2*GPR*Patm;

fprintf('Compressing fuel\n');
[PfuelComp,TfuelComp,state_2f,path_1f_2f] = compTurb(fuel,P1f,T1f,P2f,eta_fuelComp,2*nsteps);
h2f = state_2f.H;
T2f = state_2f.T;
%% set up water
fluid = Water;
P0w = Patm;
T0w = Tatm;
set(fluid,'P',P0w,'T',T0w);
h0w = enthalpy_mass(fluid);
%% compress water
P2w = P2f;
P1w = P2w/(econ_pratio*boil_pratio*superheat_pratio);
[PwaterPump,TwaterPump,state_1w,path_0w_1w] = compTurb(fluid,P0w,T0w,P1w,eta_feedPump,nsteps);
T1w = state_1w.T;
h1w = state_1w.H;
%% Mix air-fuel-steam and burn it 
%  Steam is generated through an HRSG

fprintf('Doing some rubbish with burning fuel and injecting steam..!\n');

SAratioVec = 0.0:0.025:0.4;
SAratioVec = 0.4;

TpinchVec = zeros(length(SAratioVec),1);
TurbineExVec = TpinchVec;
HRSGExVec = TpinchVec;
TsteamInjVec = TpinchVec;
Tturb_steamVec = TpinchVec;
mdot_fuelVec = TpinchVec;
TITachVec = TpinchVec;

LHVfuel = LHV_mass(fuel);
LHVeffVec = TpinchVec;
AirSpecificWork = TpinchVec;

% guess
mdot_fuel = 0.05; % solve for this iteratively

% run dry to get the initial guess for temperatures
P3m = P2a;
steam = GRI30;
% get initial guess for TIT
[T3m,T4m,h4m,x3m,x4m] = mixerburner(gasmix,fuel,air,steam,mdot_fuel,mdot_air,0,...
        P3m,burner_pratio);

% run the dry gas mixture through the turbine
P4m = pressure(gasmix);
P5m = Patm;
[PgasTurb,TgasTurb,state_5m,path_4m_5m] = compTurb(gasmix,P4m,T4m,P5m,eta_gasTurb,nsteps);
T5m = state_5m.T;
% use this dry T5m as the initial guess for the steam temperature
% further iterations use the previous T5m lol

for iter = 1:1:length(SAratioVec)
    
    SAratio = SAratioVec(iter); % parameter sweep this
    
    
    mdot_steam = mdot_air*SAratio;
    mdot_mix = mdot_air + mdot_fuel + mdot_steam;
    
    
    TIThelper = @(mdot_fuel)(TITfinder(mdot_fuel,gasmix,air,fuel,fluid,steam,T5m,P2w,T3m,P3m,x3m,...
        TIT,Patm,eta_gasTurb,nsteps,mdot_air,mdot_steam,T1w,P1w,...
        HRSG_eff,econ_pratio,boil_pratio,superheat_pratio,burner_pratio)...
        );
    
    opts = optimset('TolX',1e-7);
    mdot_fuel = fzero(TIThelper,mdot_fuel,opts); 
    
    [dTIT,retStr] = TIThelper(mdot_fuel);
    

    waterStates = retStr.waterStates;
    mixStates = retStr.mixStates;
    T3m = retStr.T3m;
    mf3m = retStr.x3m;
    T4m = retStr.T4m;
    h4m = retStr.h4m;
    mf4m = retStr.x4m;
    T5m = retStr.T5m;
    h5m = retStr.h5m;
    mf5m = retStr.x5m;

    state_5m = retStr.state_5m;
    path_4m_5m = retStr.path_4m_5m;
    TpinchActual = retStr.TpinchActual;
    TITachieved = retStr.TITachieved;

    % compute the work output from the turbine
    Wturb = mdot_mix*(h4m - h5m);
    WcompAir = mdot_air*(h2a - h1a);
    WcompFuel = mdot_fuel*(h2f - h1f);
    WfeedPump = mdot_steam*(h1w - h0w);

    Wnet = Wturb - (WcompAir + WcompFuel + WfeedPump);

    % disp(TpinchActual);
    fprintf('%3.3f\t%3.3f\t%3.3f\n',SAratio,TITachieved,TpinchActual);

    TITachVec(iter) = TITachieved;
    TpinchVec(iter) = TpinchActual;
    TurbineExVec(iter) = T5m;
    HRSGExVec(iter) = mixStates.m4.T;
    TsteamInjVec(iter) = waterStates.w4.T;
    Tturb_steamVec(iter) = T5m - waterStates.w4.T;
    mdot_fuelVec(iter) = mdot_fuel;
    LHVeffVec(iter) = Wnet/(mdot_fuel*LHVfuel);
    AirSpecificWork(iter) = Wnet/mdot_air;
end
%% make plots
if length(SAratioVec) > 2
    figure(1);clf;
    plot(SAratioVec,TurbineExVec,'b-','LineWidth',2.0);
    hold on;
    plot(SAratioVec,HRSGExVec,'-','LineWidth',2.0,'Color',[0 0.5 0]);
    plot(SAratioVec,TsteamInjVec,'r-','LineWidth',2.0);
    plot(SAratioVec,Tturb_steamVec,'-','LineWidth',2.0,'Color',[0 0.85 0.85]);
    plot(SAratioVec,TpinchVec,'m-','LineWidth',2.0);
    
    xlabel('Steam/Air Mass Ratio');
    ylabel('T or \DeltaT (K)');
    
    lgd = legend('T Turbine Exhaust','T HRSG Exhaust', ...
        'T Steam Injection','\DeltaT Turbine-Steam','DeltaT Pinch-Point');
    
    set(lgd,'Location','Best')
    set(gca,'FontSize',26,'FontWeight','bold');
    set(gcf,'Position',[680 102 1180 776],'Color','White');
    
    h = gcf;
    saveas(h,'p4_tempPlots.fig');
    print(gcf,'p4_tempPlots.png','-dpng','-r300'); 
    
    % First law efficiency plot
    figure(2);clf;
    plot(SAratioVec,100*LHVeffVec,'k-','LineWidth',2);
    xlabel('Steam/Air Mass Ratio');
    ylabel('LHV efficiency (%)');
    set(gca,'FontSize',26,'FontWeight','bold');
    set(gcf,'Position',[680 102 1180 776],'Color','White');
    
    h = gcf;
    saveas(h,'p4_LHVeff.fig');
    print(gcf,'p4_LHVeff.png','-dpng','-r300'); 
    
    % Air specific work
    figure(3);clf;
    plot(SAratioVec,AirSpecificWork/1e3,'k-','LineWidth',2);
    xlabel('Steam/Air Mass Ratio');
    ylabel('Air specific work (kJ/kg-air)');
    set(gca,'FontSize',26,'FontWeight','bold');
    set(gcf,'Position',[680 102 1180 776],'Color','White');
    
    h = gcf;
    saveas(h,'p4_airSpecWork.fig');
    print(gcf,'p4_airSpecWork.png','-dpng','-r300'); 

end
%% Exergy plots
exergyGasPlaceHolder = GRI30;
exergyWaterPlaceHolder = Water;
exergyFun = @(canteraObject)(exergy_mass(canteraObject,Tatm,Patm,xair));
exergyFlowFun = @(canteraObject)(flowExergy_mass(canteraObject,Tatm,Patm,xair));

% gas states zero --> fuel and air dead states
set(exergyGasPlaceHolder,'T',Tatm,'P',Patm,'X',xair);
x0a = exergyFun(exergyGasPlaceHolder);
xf0a = exergyFlowFun(exergyGasPlaceHolder);
set(exergyGasPlaceHolder,'T',Tatm,'P',Patm,'X',xfuel);
x0f = exergyFun(exergyGasPlaceHolder);
xf0f = exergyFlowFun(exergyGasPlaceHolder);

% gas states one --> fuel and air
set(exergyGasPlaceHolder,'T',T1a,'P',P1a,'X',xair);
x1a = exergyFun(exergyGasPlaceHolder);
xf1a = exergyFlowFun(exergyGasPlaceHolder);
set(exergyGasPlaceHolder,'T',T1a,'P',P1a,'X',xfuel);
x1f = exergyFun(exergyGasPlaceHolder);
xf1f = exergyFlowFun(exergyGasPlaceHolder);

% gas states two --> fuel and air
set(exergyGasPlaceHolder,'T',T2a,'P',P2a,'X',xair);
x2a = exergyFun(exergyGasPlaceHolder);
xf2a = exergyFlowFun(exergyGasPlaceHolder);
set(exergyGasPlaceHolder,'T',T2f,'P',P2f,'X',xfuel);
x2f = exergyFun(exergyGasPlaceHolder);
xf2f = exergyFlowFun(exergyGasPlaceHolder);

% gas states three, four, five and six --> only mix
set(exergyGasPlaceHolder,'T',T3m,'P',P3m,'X',mf3m);
x3m = exergyFun(exergyGasPlaceHolder);
xf3m = exergyFlowFun(exergyGasPlaceHolder);

set(exergyGasPlaceHolder,'T',T4m,'P',P4m,'X',mf4m);
x4m = exergyFun(exergyGasPlaceHolder);
xf4m = exergyFlowFun(exergyGasPlaceHolder);

set(exergyGasPlaceHolder,'T',T5m,'P',P5m,'X',mf5m);
x5m = exergyFun(exergyGasPlaceHolder);
xf5m = exergyFlowFun(exergyGasPlaceHolder);

set(exergyGasPlaceHolder,'T',mixStates.m4.T,'P',mixStates.m4.P,'X',mixStates.m4.x);
x6m = exergyFun(exergyGasPlaceHolder);
xf6m = exergyFlowFun(exergyGasPlaceHolder);




% water state zero
set(exergyWaterPlaceHolder,'T',Tatm,'P',Patm);
x0w = exergyFun(exergyWaterPlaceHolder);
xf0w = exergyFlowFun(exergyWaterPlaceHolder);

% water state one
set(exergyWaterPlaceHolder,'T',T1w,'P',P1w);
x1w = exergyFun(exergyWaterPlaceHolder);
xf1w = exergyFlowFun(exergyWaterPlaceHolder);

% water state two
set(exergyWaterPlaceHolder,'T',waterStates.w4.T,'P',waterStates.w4.P);
x2w = exergyFun(exergyWaterPlaceHolder);
xf2w = exergyFlowFun(exergyWaterPlaceHolder);

% water state three
set(exergyGasPlaceHolder,'T',T3m,'P',waterStates.w4.P,'X','H2O:1');
x3w = exergyFun(exergyGasPlaceHolder);
xf3w = exergyFlowFun(exergyGasPlaceHolder);

mdot_w = mdot_steam;
mdot_mix = mdot_air + mdot_fuel + mdot_steam;
W_fuel_compressor = WcompFuel;
W_air_compressor = WcompAir;
W_gross_gas_turbine = Wturb;
W_feed_pump = WfeedPump;
Eff_LHV = LHVeffVec(end);

Eff_exergy = nan;

%% gas mixer and burner function
function [T3m,T4m,h4m,x3m,x4m] = mixerburner(gasmix,fuel,air,steam,mdot_fuel,mdot_air,mdot_steam,...
    P3m,burner_pratio)
    mdot_mix = mdot_fuel + mdot_air + mdot_steam;
    
    xfuel = moleFractions(fuel);
    xair = moleFractions(air);
    xsteam = moleFractions(steam);
    
    MMWfuel = meanMolecularWeight(fuel);
    MMWair = meanMolecularWeight(air);
    MMWsteam = meanMolecularWeight(steam);

    h3m = (mdot_fuel*enthalpy_mass(fuel) + mdot_air*enthalpy_mass(air) ...
        + mdot_steam*enthalpy_mass(steam))/mdot_mix;
    x3m = (mdot_fuel/MMWfuel)*xfuel + (mdot_air/MMWair)*xair + (mdot_steam/MMWsteam)*xsteam;
    x3m = x3m/sum(x3m);

    % set the gas mixture
    set(gasmix,'H',h3m,'P',P3m,'X',x3m);
    T3m = temperature(gasmix);
    
    % Burn the gas mixture
    P4m=burner_pratio*P3m;
    set(gasmix,'H',h3m,'P',P4m,'X',x3m);
    equilibrate(gasmix,'HP');
    T4m = temperature(gasmix);
    h4m = enthalpy_mass(gasmix);
    x4m = moleFractions(gasmix);
end

%% iteratively solve this for TIT 
function [dTIT,retStr] = TITfinder(mdot_fuel,gasmix,air,fuel,fluid,steam,T5m,P2w,T3m,P3m,x3m,...
    TIT,Patm,eta_gasTurb,nsteps,mdot_air,mdot_steam,T1w,P1w,...
    HRSG_eff,econ_pratio,boil_pratio,superheat_pratio,burner_pratio)
    
    set(steam,'T',T5m,'P',P2w,'X','H2O:1');
    prevTemp = 0;
    currTemp = T5m;
    
    % these iterations are required because 
    % steam injection temp affects the turbine exhaust temp
    % which again impacts the steam temp through the HRSG
    while abs(prevTemp - currTemp)>1e-1
        % burn gases
        prevTemp = currTemp;
        set(gasmix,'T',T3m,'P',P3m,'X',x3m);
        [T3m,T4m,h4m,x3m,x4m] = mixerburner(gasmix,fuel,air,steam,mdot_fuel,mdot_air,mdot_steam,...
            P3m,burner_pratio);
        TITactual = T4m;
    
        % run the gas mixture through the turbine
        P4m = pressure(gasmix);
        P5m = Patm;
        [PgasTurb,TgasTurb,state_5m,path_4m_5m] = compTurb(gasmix,P4m,T4m,P5m,eta_gasTurb,nsteps);
        T5m = state_5m.T;
        equilibrate(gasmix,'TP');
        x5m = moleFractions(gasmix);
        h5m = state_5m.H;

        % run through HRSG
        set(fluid,'T',T1w,'P',P1w);
        set(gasmix,'T',T5m,'P',P5m,'X',x5m);
        mdot_mix = mdot_air+mdot_fuel+mdot_steam;
        [Tpinch,waterStates,mixStates] = HRSG(gasmix,fluid,mdot_mix,mdot_steam,...
            HRSG_eff,econ_pratio,boil_pratio,superheat_pratio);
        set(steam,'T',waterStates.w4.T,'P',waterStates.w4.P);
        currTemp = waterStates.w4.T;
        % disp(currTemp);
        TpinchActual = mixStates.m3.T - waterStates.w2.T;
    end
    dTIT = TITactual - TIT;

    retStr.waterStates = waterStates;
    retStr.mixStates = mixStates;
    retStr.T3m = T3m;
    retStr.x3m = x3m;
    retStr.T4m = T4m;
    retStr.h4m = h4m;
    retStr.x4m = x4m;
    retStr.T5m = T5m;
    retStr.x5m = x5m;
    retStr.h5m = h5m;
    retStr.state_5m = state_5m;
    retStr.path_4m_5m = path_4m_5m;
    retStr.TpinchActual = TpinchActual;
    retStr.TITachieved = TITactual;
end