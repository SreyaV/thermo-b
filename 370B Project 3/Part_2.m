% Working file for simple-cycle steam turbine engine.
% C.F. Edwards, 1/9/10

clear all
format compact
fprintf('\n****************************************************************\n')

% Set step size for compression/expansion.
steps = 500;

% Cycle specifications:
To = 293.15         % K
Po = 100000         % Pa
Pstorage = 10e5     % Pa
Pcond = 6800        % Pa
Economizer_Pressure_Ratio             = 0.92
Evaporator_Pressure_Ratio                 = 0.92
Superheater_Pressure_Ratio            = 0.92
Condenser_Pressure_Ratio              = 0.92
Steam_Turbine_Polytropic_Efficiency   = 0.80
Condensate_Pump_Polytropic_Efficiency = 0.85
Feed_Pump_Polytropic_Efficiency       = 0.85
Turbine_Exit_Quality                  = 0.85
Deaerator_Pressure_Ratio              = 0.92
Feed_Pump_Pressure_Ratio              = (170e5/Pstorage)/Economizer_Pressure_Ratio

% Use the following lines to remove the feedwater-heating effect of the
% deaerator.  Comment them out for normal analysis.
% Pstorage = Pcond;
% Feed_Pump_Pressure_Ratio = (40e5/Pcond)/Economizer_Pressure_Ratio;
% Deaerator_Pressure_Ratio = 1

% Define a mass flow rate so that extensive values are available if
% desired.  (Set to unity if not of interest.)
mdot_water = 1;

% Process some water for a Rankine bottoming cycle
water = Solution('liquidvapor.yaml','water');
Pc = critPressure(water);
Tc = critTemperature(water);
vc = 1/critDensity(water);
set(water,'T',Tc,'P',Pc);
sc = entropy_mass(water);
Tt = 273.16;

% Make a vapor dome
dT = (Tc-0.1-Tt)/100;
i = 1;
for T=Tt:dT:Tc
    Tsatline(i) = T;
    setState_Tsat(water,[T 0]);
    sliqline(i) = entropy_mass(water);
    setState_Tsat(water,[T 1]);
    svapline(i) = entropy_mass(water);
    i = i+1;
end
% Add the critical point at the top.
Tsatline(i) = Tc;
sliqline(i) = sc;
svapline(i) = sc;
    
% Start with saturated liquid water at storage conditions.
Pw1 = Pstorage;
setState_Psat(water,[Pw1 0]);
Tw1 = temperature(water);
sw1 = entropy_mass(water);
hw1 = enthalpy_mass(water);

% The second water point is at high pressure after the feed pump.
% Walk up in pressure adjusting via the polytropic efficiency.
Pw2 = Feed_Pump_Pressure_Ratio*Pw1;
dP  = (Pw2 - Pw1)/steps;
s   = sw1;
h   = hw1;
for P = Pw1:dP:Pw2
    % Find the isentropic state for the new pressure.
    set(water,'S',s,'P',P);
    % Get the isentropic enthalpy difference.
    hs = enthalpy_mass(water);
    dhs = hs - h; 
    % The actual difference is larger due to inefficiency.
    h = h + dhs/Feed_Pump_Polytropic_Efficiency;
    % Find the new state.
    set(water,'H',h,'P',P);
    % Get the entropy for the next step in the process.
    s = entropy_mass(water);
end
Tw2 = temperature(water);
hw2 = h;
sw2 = s;

% Find the isentropic efficiency for this pressure ratio just for fun.
set(water,'S',sw1,'P',Pw2);
hw2s = enthalpy_mass(water);
Feed_Pump_Isentropic_Efficiency = (hw2s - hw1)/(hw2 - hw1)

% Do the economizer.
Pw3 = Pw2*Economizer_Pressure_Ratio;
setState_Psat(water,[Pw3 0]);
Tw3 = temperature(water);
sw3 = entropy_mass(water);
hw3 = enthalpy_mass(water);

% Go across the steam drum.
Pw4 = Pw3*Evaporator_Pressure_Ratio;
setState_Psat(water,[Pw4 1]);
Tw4 = temperature(water);
sw4 = entropy_mass(water);
hw4 = enthalpy_mass(water);

% Find the peak temperature required to meet an exit quality spec.
% The turbine outlet state is known, and the inlet pressure is known.
Pw5 = Pw4*Superheater_Pressure_Ratio
Pw6 = Pcond
setState_Psat(water,[Pw6 Turbine_Exit_Quality]);
Tw6 = temperature(water);
sw6 = entropy_mass(water);
hw6 = enthalpy_mass(water);
vw6 = 1/density(water);

% There are a couple of ways to do the expansion.  The most obvious is to
% step down in pressure.  By marching through a series of starting
% temperatures, you can find the exit quality.  Another choice is to
% march in specific volume (or density), solving for the pressure
% accompanying each step while still iterating the temperature as an outer
% variable.  The last choice is the least obvious: walk back up the
% polytrope from the known State 6 back up to State 5.  More about this one
% below.

method = 1;
switch method
    case 1
        % Pressure step sizes are tricky because Cantera has trouble finding the right
        % density as you pass the critical temperature.  Use try/catch/end to get past
        % the bad points.  Usually full convergence takes ~500 steps, .5 K.
        dP  = (Pw5 - Pw6)/steps;
        dT  = 1;
        % Start from State 4 and walk upwards.
        Tstart = Tw4;
        for Tw5=Tstart:dT:1000
            Tw5;
            try
            % Get the starting state for this temperature.
            set(water,'T',Tw5,'P',Pw5);
            sw5 = entropy_mass(water);
            hw5 = enthalpy_mass(water);
            vw5 = 1/density(water);
            % Expand the steam to the condenser pressure.
            % Walk down in pressure adjusting via the polytropic efficiency.
            s = sw5;
            h = hw5;
            T = Tw5;
            v = vw5;
            for P=Pw5:-dP:Pw6
                % Find the isentropic state for the new pressure.
                set(water,'S',s,'P',P);
                % Get the isentropic enthalpy difference.
                hs = enthalpy_mass(water);
                dhs = h - hs;
                % The actual difference is smaller due to inefficiency.
                h = h - dhs*Steam_Turbine_Polytropic_Efficiency;
                % Find the new state.
                set(water,'H',h,'P',P);
                % Get the entropy for the next step in the process.
                s = entropy_mass(water);
                % Get the temperature.  Print this to check for near-critical
                % problems.
                T = temperature(water);
                v = 1/density(water);
            end
            Tw6 = temperature(water);
            hw6 = h;
            sw6 = s;
            % Check the quality to see if you crossed yet.
            Quality = vaporFraction(water)
            if(Quality >= Turbine_Exit_Quality)
                break
            end
            catch
                disp('Jumping critical isotherm...')
                water = importPhase('liquidvapor.xml','water');
            end
        end
        Tw5
        Quality
        
    case 2
        % Solve by specific volume.  Steps can be tricky here too.
        % Find the specific volume at the exit state.
        dT  = 1;
        % Start from State 4 and walk upwards.
        Tstart = Tw4;
%         Tstart = 790;
        for Tw5=Tstart:dT:1000
            % Get the starting state for this temperature.
            set(water,'T',Tw5,'P',Pw5);
            sw5 = entropy_mass(water);
            hw5 = enthalpy_mass(water);
            vw5 = 1/density(water);
            dv = (vw6-vw5)/steps;
            % Expand the steam to the condenser state.
            % Walk in specific volume adjusting via the polytropic efficiency.
            s = sw5;
            h = hw5;
            T = Tw5;
            for v = vw5:dv:vw6
                % Find the isentropic state for the new pressure.  This state
                % can still run into trouble by hitting critical T.
                try
                    set(water,'S',s,'V',v);
                    % Get the isentropic enthalpy difference and pressure.
                    P = pressure(water);
                    hs = enthalpy_mass(water);
                    dhs = h - hs;
                    % The actual difference is smaller due to inefficiency.
                    h = h - dhs*Steam_Turbine_Polytropic_Efficiency;
                    % Find the new state.  It must have the same pressure or
                    % you will violate the definition of polytropic efficiency.
                    set(water,'H',h,'P',P);
                    % Get the entropy for the next step in the process.
                    s = entropy_mass(water);
                    % Get the temperature.  Print this to check for near-critical
                    % problems.
                    T = temperature(water);
                catch
                    disp('Jumping critical isotherm...')
                    water = importPhase('liquidvapor.xml','water');
                end
            end
            Tw6 = temperature(water);
            hw6 = h;
            sw6 = s;
            % Check the quality to see if you crossed yet.
            Quality = vaporFraction(water)
            if(Quality >= Turbine_Exit_Quality)
                break
            end
        end
        Tw5
        Quality
       
    case 3
        % Do pressure stepping but come from below.  Essentially you are
        % doing the problem backwards. Need to be fully converged to get
        % the correct answer.  (I suppose that is obvious.)  The point here
        % is that if you do this in discreet steps, the definition of
        % polytropic efficiency appears to be violated.  But when the steps
        % become small enough, life becomes linear, and you get the correct
        % path to match the definition.  (And the problem is explicit.)
        dP  = (Pw5 - Pw6)/steps;
        
        % Reverse-expand the steam to the turbine inlet pressure.
        % Walk up in pressure adjusting via the polytropic efficiency.
        s = sw6;
        h = hw6;
        T = Tw6;
        for P=Pw6:dP:Pw5
            % Find the isentropic state for the new pressure.
            try
                set(water,'S',s,'P',P);
                % Get the isentropic enthalpy difference.
                hs = enthalpy_mass(water);
                dhs = h - hs;
                % The actual difference is smaller due to inefficiency.
                h = h - dhs*Steam_Turbine_Polytropic_Efficiency;
                % Find the new state.
                set(water,'H',h,'P',P);
                % Get the entropy for the next step in the process.
                s = entropy_mass(water);
                % Get the temperature.  Print this to check for near-critical
                % problems.
                T = temperature(water);
            catch
                disp('Jumping critical isotherm...')
                water = importPhase('liquidvapor.xml','water');
            end
        end
        Tw5 = T
        sw5 = s;
        hw5 = h;
end

% Find the isentropic efficiency for this pressure ratio just for fun.
set(water,'S',sw5,'P',Pw6);
hw6s = enthalpy_mass(water);
Steam_Turbine_Isentropic_Efficiency = (hw5 - hw6)/(hw5 - hw6s)

% Condense the steam back to saturated water.
Pw7 = Pcond*Condenser_Pressure_Ratio;
setState_Psat(water,[Pw7 0]);
Tw7 = temperature(water);
sw7 = entropy_mass(water);
hw7 = enthalpy_mass(water);

% The last water point is before the deaerator via a condensate pump.
% Walk up in pressure adjusting via the polytropic efficiency.
Pw8 = Pw1/Deaerator_Pressure_Ratio;
dP  = (Pw8 - Pw7)/steps;
s   = sw7;
h   = hw7;
for P = Pw7:dP:Pw8
    set(water,'S',s,'P',P);
    hs = enthalpy_mass(water);
    dhs = hs - h;
    h = dhs/Condensate_Pump_Polytropic_Efficiency + h;
    set(water,'H',h,'P',P);
    s = entropy_mass(water);
end
Tw8 = temperature(water);
hw8 = h;
sw8 = s;

% Find the isentropic efficiency for this pressure ratio just for fun.
set(water,'S',sw7,'P',Pw8);
hw8s = enthalpy_mass(water);
Condensate_Pump_Isentropic_Efficiency = (hw8s - hw7)/(hw8 - hw7)

% Find the state at the pressure of the extraction.
% Note that both the liquid from the condenser and the steam extracted have
% the same inlet pressure Pw8, and same outlet pressure Pw1, for the
% dearator.  (So the pressure drop is "across the box" from ins to out.)
Pwe = Pw8;

% Expand the steam from State 5 to the extraction pressure.  We don't yet
% know how much we need to extract, we just have to find the state.
% Walk down in pressure adjusting via the polytropic efficiency.
s   = sw5;
h   = hw5;
dP  = (Pw5 - Pwe)/steps;
missed_on_last = false;
for P = Pw5:-dP:Pwe
    try
        set(water,'S',s,'P',P);
        hs = enthalpy_mass(water);
        dhs = h - hs;
        h = h - dhs*Steam_Turbine_Polytropic_Efficiency;
        set(water,'H',h,'P',P);
        s = entropy_mass(water);
        missed_on_last = false;
    catch
        disp('Jumping critical isotherm...')
        water = importPhase('liquidvapor.xml','water');
    end
end
Twe = temperature(water);
hwe = h;
swe = s;

% Find the fraction of water that must be extracted for the deaerator.
if Pw1==Pwe
    % If the extract pressure is the same as the deaerator outlet, don't
    % extract any.
    f = 0
else
    % Otherwise use an energy balance on the deaerator to find the amount:
    % mdot1*h1 = mdote*he + mdot8*h8
    % h1 = (mdote/mdot1)*he + (mdot8/mdot1)*h8  Let f = mdote/mdot1
    % h1 = f*he + (1-f)*h8
    % h1-h8 = f*(he-h8)
    % f = (h1-h8)/(he-h8)
    f = (hw1-hw8)/(hwe-hw8)
end

% Make an array of plotting points.
Tpts = [Tw1 Tw2 Tw3 Tw4 Tw5 Twe Tw6 Tw7 Tw8 Tw1];
spts = [sw1 sw2 sw3 sw4 sw5 swe sw6 sw7 sw8 sw1];

% Make an array that shows the extraction.  
% Must find the saturated vapor point.
setState_Psat(water,[Pwe 1]);
Twesat1 = temperature(water);
swesat1 = entropy_mass(water);
setState_Psat(water,[Pwe 0]);
Twesat0 = temperature(water);
swesat0 = entropy_mass(water);
% We don't know the path for the extracted steam.  It might be better to
% show it as a straight line to the deaerator state.  But in the lines
% below, the saturation point is included (as though the pressure does not
% drop until after you cross the dome, which is not true).  The reason for
% depicting it this way is that it looks very odd when the extraction is up
% in the superheat region and you get a big diagonal on the diagram.
% Tptse = [Twe Twesat1 Twesat0 Tw1];
% sptse = [swe swesat1 swesat0 sw1];
% Here it is the other way.  Comment one of these out.
Tptse = [Twe Tw1];
sptse = [swe sw1];

% Find the various energy transfers.
W_pump_condensate     = mdot_water*(1-f)*(hw8 - hw7)
W_pump_feedwater      = mdot_water*(hw2 - hw1)
Q_economizer          = mdot_water*(hw3 - hw2)
Q_evaporator              = mdot_water*(hw4 - hw3)
Q_superheater         = mdot_water*(hw5 - hw4)
W_gross_steam_turbine = mdot_water*(hw5 - hwe) + mdot_water*(1-f)*(hwe - hw6)
W_net_steam_turbine   = W_gross_steam_turbine - W_pump_condensate - W_pump_feedwater
Q_condenser           = mdot_water*(1-f)*(hw6 - hw7)
Eff_steam_turbine     = W_net_steam_turbine/(Q_economizer + Q_evaporator + Q_superheater)
Specific_Work         = W_net_steam_turbine/mdot_water

% Make T-s plot.
figure(1)
clf
hold on
plot(sliqline/1000,Tsatline,'k')
plot(svapline/1000,Tsatline,'k')
plot(spts/1000,Tpts,'bo--')
plot(sptse/1000,Tptse,'b--')
text(sw1/1000,Tw1,'1')
text(sw2/1000,Tw2,'2')
text(sw3/1000,Tw3,'3')
text(sw4/1000,Tw4,'4')
text(sw5/1000,Tw5,'5')
text(sw6/1000,Tw6,'6')
text(sw7/1000,Tw7,'7')
text(sw8/1000,Tw8,'8')
text(swe/1000,Twe,'e')
text(2.75,820,sprintf('Thermal Efficiency:   %.1f%%',100*Eff_steam_turbine))
text(2.75,760,sprintf('Specific Work:         %.2f MJ/kg',Specific_Work/1e6))
text(2.75,700,sprintf('Max Temperature:     %.0f K',Tw5))
xlabel('Specific Entropy (kJ/kg-K)')
ylabel('Temperature (K)')
hold off
scale = axis;
axis([2 scale(2) scale(3) 900])
plotfixer
legend('off')
