% Make a T-s, u-s, and x-s diagrams for water.
% C.F. Edwards, 1/8/11

clear all
format compact
fprintf('\n********************************************************\n')

% Set the environmental state
global To Po mu_o

% Set the environmental state.
% Make a Cantera gas object and get indices for species of interest.
gas  = Solution('gri30.yaml','gri30');
N    = nSpecies(gas);
iO2  = speciesIndex(gas,'O2');
iN2  = speciesIndex(gas,'N2');
iCO2 = speciesIndex(gas,'CO2');
iH2O = speciesIndex(gas,'H2O');
iAR  = speciesIndex(gas,'AR');
To = 25+273.15
Po = 101325
xo = zeros(1,N);
xo(iN2)  = 0.753611;
xo(iO2)  = 0.202157;
xo(iH2O) = 0.034820;
xo(iAR)  = 0.009015;
xo(iCO2) = 0.000397;
% Check for sum to unity.
Sum_xo = sum(xo)
if Sum_xo ~= 1
    disp('Mole fractions dont sum to unity...')
end
% Get the chemical potentials.  These will be used by the exergy function
% and are accessed through global storage.
set(gas,'T',To,'P',Po,'X',xo);
mu_o = chemPotentials(gas);

% Set up a water object to work with in Cantera.
water = Solution('liquidvapor.yaml','water');

% Get the critical point props.
Tc = critTemperature(water)
Pc = critPressure(water)
set(water,'T',Tc,'P',Pc);
sc = entropy_mass(water);
uc = intEnergy_mass(water);
xc = exergy_mass(water);

% Set the triple point props.  Use a value epsilon above Tt to get the
% bottom edge of the vapor dome.
Tt = 273.17
Pt = satPressure(water,Tt)
setState_Tsat(water,[Tt 0]);
sft = entropy_mass(water);
uft = intEnergy_mass(water);
vft = 1/density(water);
xft = exergy_mass(water);
setState_Tsat(water,[Tt 1]);
sgt = entropy_mass(water);
ugt = intEnergy_mass(water);
vgt = 1/density(water);
xgt = exergy_mass(water);

% Set the limits for data curves.
Tmin = 274;
Tmax = 1000;
Pmin = Pt;
Pmax = 1000*100000;

% Make a vapor dome.
dT = (Tc-Tmin)/100;
i = 1;
%setState_Tsat(water,[Tmin 0]);
for T=Tmin:dT:Tc-dT    % Stop short of the critical point.
    Tsatline(i) = T;
    setState_Tsat(water,[T 0]);
    sliqline(i) = entropy_mass(water);
    uliqline(i) = intEnergy_mass(water);
    xliqline(i) = exergy_mass(water);
    setState_Tsat(water,[T 1]);
    svapline(i) = entropy_mass(water);
    uvapline(i) = intEnergy_mass(water);
    xvapline(i) = exergy_mass(water);
    i = i+1;
end
Tsatline(i) = Tc;   % Add the critical point now.
sliqline(i) = sc;
uliqline(i) = uc;
xliqline(i) = xc;
svapline(i) = sc;
uvapline(i) = uc;
xvapline(i) = xc;

% Start a set of isobaric curves.
Plist = [Pmin/100000 0.1 1 10 100];  % Pressure in bar
for j=1:1:length(Plist)
    water = Solution('liquidvapor.yaml','water');
    P = Plist(j)*100000

    % Do the compressed liquid side.
    setState_Psat(water,[P 0]);
    Tsat = temperature(water);
    dT = (Tsat - Tmin)/50;
    i = 1;
    for T=Tmin:dT:Tsat-dT  % Stop short of saturation.
        Tpresline(j,i) = T;
        set(water,'T',T,'P',P);
        spresline(j,i) = entropy_mass(water);
        upresline(j,i) = intEnergy_mass(water);
        xpresline(j,i) = exergy_mass(water);
        i = i+1;
    end
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    setState_Tsat(water,[Tsat 0]);
    spresline(j,i) = entropy_mass(water);
    upresline(j,i) = intEnergy_mass(water);
    xpresline(j,i) = exergy_mass(water);
    i = i+1;

    % Now go across the dome.
    dq = 1/50;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Tpresline(j,i) = Tsat;
        setState_Psat(water,[P q]);
        spresline(j,i) = entropy_mass(water);
        upresline(j,i) = intEnergy_mass(water);
        xpresline(j,i) = exergy_mass(water);
        i = i+1;
    end
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    setState_Psat(water,[P 1]);
    spresline(j,i) = entropy_mass(water);
    upresline(j,i) = intEnergy_mass(water);
    xpresline(j,i) = exergy_mass(water);
    i = i+1;

    % Do the vapor side.
    dT = (Tmax - Tsat)/50;
    for T=Tsat+dT:dT:Tmax  % Start just above saturation.
        Tpresline(j,i) = T;
        set(water,'T',T,'P',P);
        spresline(j,i) = entropy_mass(water);
        upresline(j,i) = intEnergy_mass(water);
        xpresline(j,i) = exergy_mass(water);
        i = i+1;
    end
end

% Add an isobar above the critical pressure.
P = 1000;           % In bar
Plist = [Plist P];  % Add to the list
P = P*100000;       % In Pascals
j = length(Plist);
dT = (Tmax - Tmin)/150;
i = 1;
for T=Tmin:dT:Tmax  % Stop short of saturation.
    Tpresline(j,i) = T;
    set(water,'T',T,'P',P);
    spresline(j,i) = entropy_mass(water);
    upresline(j,i) = intEnergy_mass(water);
    xpresline(j,i) = exergy_mass(water);
    i = i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isotherms.
Tlist = [Tmin 400 500 600 700 800 900 1000];
for j=1:1:3    % Do the part of list below the critical point. 
    water = Solution('liquidvapor.yaml','water');
    T = Tlist(j)

    % Do the compressed liquid side.
    setState_Tsat(water,[T 0]);
    Psat = pressure(water);
    logdP = (log(Pmax) - log(Psat))/50;
    i = 1;
    for logP=log(Pmax):-logdP:log(Psat)+logdP  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(water,'T',T,'P',P);
        stempline(j,i) = entropy_mass(water);
        utempline(j,i) = intEnergy_mass(water);
        xtempline(j,i) = exergy_mass(water);
        i = i+1;
    end
    Ttempline(j,i) = T;   % Add the saturation point now.
    % set(water,'T',T,'P',Psat);
    setState_Tsat(water,[T 0]);
    stempline(j,i) = entropy_mass(water);
    utempline(j,i) = intEnergy_mass(water);
    xtempline(j,i) = exergy_mass(water);
    i = i+1;

    % Now go across the dome.
    dq = 1/50;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Ttempline(j,i) = T;
        setState_Psat(water,[Psat q]);
        stempline(j,i) = entropy_mass(water);
        utempline(j,i) = intEnergy_mass(water);
        xtempline(j,i) = exergy_mass(water);
        i = i+1;
    end
    Ttempline(j,i) = T; % Add the saturation point now.
    setState_Psat(water,[Psat 1]);
    stempline(j,i) = entropy_mass(water);
    utempline(j,i) = intEnergy_mass(water);
    xtempline(j,i) = exergy_mass(water);
    i = i+1;

    % Do the vapor side.
    logdP = (log(Psat) - log(Pmin))/50;
    for logP=log(Psat)-logdP:-logdP:log(Pmin)  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(water,'T',T,'P',P);
        stempline(j,i) = entropy_mass(water);
        utempline(j,i) = intEnergy_mass(water);
        xtempline(j,i) = exergy_mass(water);
        i = i+1;
    end
end

% Add isotherms above the critical temperature.
for j=j+1:1:length(Tlist)
    water = Solution('liquidvapor.yaml','water');
    T = Tlist(j)

    logdP = (log(Pmax) - log(Pmin))/150;
    i = 1;
    for logP=log(Pmax):-logdP:log(Pmin)  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(water,'T',T,'P',P);
        stempline(j,i) = entropy_mass(water);
        utempline(j,i) = intEnergy_mass(water);
        xtempline(j,i) = exergy_mass(water);
        i = i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isochores.
set(water,'T',Tmin,'P',Pmin);
vmin = 1/density(water);
vmax = vgt;
set(water,'T',Tmax,'V',vmin);
vlist = [0.001 0.01 0.1 1 10 100];
for j=1:1:length(vlist)
    v = vlist(j)
    dT = (Tmax - Tmin)/150;
    i = 1;
    for T=Tmin:dT:Tmax  % Stop short of saturation.
        Tvoluline(j,i) = T;
        set(water,'T',T,'V',v);
        svoluline(j,i) = entropy_mass(water);
        uvoluline(j,i) = intEnergy_mass(water);
        xvoluline(j,i) = exergy_mass(water);
        i = i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isoergs.
ulist = [-13 -13.25 -13.5 -14 -14.5 -15 -15.5]*1e6;
vmax  = [100    100   100 100   100  70    40];
vmin = .001001;
for j=1:1:length(ulist)
    u = ulist(j)
    lvhigh = log(vmax(j));
    lvlow  = log(vmin);
    ldv = (lvhigh - lvlow)/300;
    i = 1;
    for lv=lvhigh:-ldv:lvlow  % Stop short of saturation.
        v = exp(lv);
        try
            setState_UV(water,[u,v]);
            Terguline(j,i) = temperature(water);
            serguline(j,i) = entropy_mass(water);
            xerguline(j,i) = exergy_mass(water);
            if(Terguline(j,i) > Tmax)
                Terguline(j,i) = Terguline(j,i-1);
                serguline(j,i) = serguline(j,i-1);
                xerguline(j,i) = xerguline(j,i-1);
            end
        catch
            disp('Trouble finding isoerg...')
            water = Solution('liquidvapor.yaml','water');
            Terguline(j,i) = Terguline(j,i-1);
            serguline(j,i) = serguline(j,i-1);
            xerguline(j,i) = xerguline(j,i-1);
        end
        i = i+1;
    end
end

% Make the plots.
figure(1)
clf
hold on
i = 1;
plot(spresline(i,:)/1000,Tpresline(i,:),'b')
plot(svoluline(i,:)/1000,Tvoluline(i,:),'-','Color',[0 .6 0])
plot(serguline(i,:)/1000,Terguline(i,:),'r')
plot(sc/1000,Tc,'kd')
plot([sft/1000 sgt/1000],[Tt Tt],'ko--')
for i=1:1:length(Plist)
    plot(spresline(i,:)/1000,Tpresline(i,:),'b')
end
for i=1:1:length(vlist)
    plot(svoluline(i,:)/1000,Tvoluline(i,:),'-','Color',[0 .6 0])
end
for i=1:1:length(ulist)
    plot(serguline(i,:)/1000,Terguline(i,:),'r')
end
plot(sliqline/1000,Tsatline,'k')
plot(svapline/1000,Tsatline,'k')
plot([sft/1000 sgt/1000],[Tt Tt],'ko--')
hold off
xlabel('Specific Entropy (kJ/kg-K)')
ylabel('Temperature (K)')
for i=1:1:length(Plist)
    text(spresline(i,150)/1000,Tpresline(i,150),num2str(Plist(i)),'Color','b')
end
for i=1:1:length(vlist)
    text(svoluline(i,150)/1000,Tvoluline(i,150),num2str(vlist(i)),'Color',[0 .6 0])
end
for i=1:1:length(ulist)
    text(serguline(i,1)/1000,Terguline(i,1),num2str(ulist(i)/1e6),'Color','r')
end
legend('Isobar (bar)','Isochore (m^3/kg)','IsoIntE (MJ/kg)','Critical Point','Triple Line')
title('Water')
plotfixer

figure(2)
clf
hold on
i = 1;
plot(spresline(i,:)/1000,upresline(i,:)/1e6,'b')
plot(stempline(i,:)/1000,utempline(i,:)/1e6,'r')
plot(svoluline(i,:)/1000,uvoluline(i,:)/1e6,'g')
plot(sc/1000,uc/1e6,'kd')
plot([sft/1000 sgt/1000],[uft/1e6 ugt/1e6],'ko--')
for i=1:1:length(Plist)
    plot(spresline(i,:)/1000,upresline(i,:)/1e6,'b')
end
for i=1:1:length(Tlist)
    plot(stempline(i,:)/1000,utempline(i,:)/1e6,'r')
end
for i=1:1:length(vlist)
    plot(svoluline(i,:)/1000,uvoluline(i,:)/1e6,'g')
end
plot(sliqline/1000,uliqline/1e6,'k')
plot(svapline/1000,uvapline/1e6,'k')
plot([sft/1000 sgt/1000],[uft/1e6 ugt/1e6],'ko--')
hold off
xlabel('Specific Entropy (kJ/kg-K)')
ylabel('Specific Internal Energy (MJ/kg)')
for i=1:1:length(Plist)
    text(spresline(i,150)/1000,upresline(i,150)/1e6,...
        num2str(Plist(i)),'Color','b')
end
for i=1:1:length(Tlist)
    text(stempline(i,150)/1000,utempline(i,150)/1e6,...
        num2str(Tlist(i)),'Color','r')
end
for i=1:1:length(vlist)
    text(svoluline(i,150)/1000,uvoluline(i,150)/1e6,...
        num2str(vlist(i)),'Color',[0 .6 0])
end
legend('Isobar (bar)','Isotherm (K)','Isochore (m^3/kg)',...
    'Critical Point','Triple Line')
title('Water')
plotfixer

figure(3)
clf
hold on
i = 1;
plot(spresline(i,:)/1000,xpresline(i,:)/1e6,'b')
plot(stempline(i,:)/1000,xtempline(i,:)/1e6,'r')
plot(svoluline(i,:)/1000,xvoluline(i,:)/1e6,'g')
plot(serguline(i,:)/1000,xerguline(i,:)/1e6,'m')
plot(sc/1000,xc/1e6,'kd')
plot([sft/1000 sgt/1000],[xft/1e6 xgt/1e6],'ko--')
for i=1:1:length(Plist)
    plot(spresline(i,:)/1000,xpresline(i,:)/1e6,'b')
end
for i=1:1:length(Tlist)
    plot(stempline(i,:)/1000,xtempline(i,:)/1e6,'r')
end
for i=1:1:length(vlist)
    plot(svoluline(i,:)/1000,xvoluline(i,:)/1e6,'g')
end
for i=1:1:length(ulist)
    plot(serguline(i,:)/1000,xerguline(i,:)/1e6,'m')
end
plot(sliqline/1000,xliqline/1e6,'k')
plot(svapline/1000,xvapline/1e6,'k')
plot([sft/1000 sgt/1000],[xft/1e6 xgt/1e6],'ko--')
hold off
xlabel('Specific Entropy (kJ/kg-K)')
ylabel('Specific Internal Exergy (MJ/kg)')
for i=1:1:length(Plist)
    text(spresline(i,150)/1000,xpresline(i,150)/1e6,...
        num2str(Plist(i)),'Color','b')
end
for i=1:1:length(Tlist)
    text(stempline(i,150)/1000,xtempline(i,150)/1e6,...
        num2str(Tlist(i)),'Color','r')
end
for i=1:1:length(vlist)
    text(svoluline(i,150)/1000,xvoluline(i,150)/1e6,...
        num2str(vlist(i)),'Color',[0 .6 0])
end
for i=1:1:length(ulist)
    text(serguline(i,1)/1000,xerguline(i,1)/1e6,...
        num2str(ulist(i)/1e6),'Color','m')
end
legend('Isobar (bar)','Isotherm (K)','Isochore (m^3/kg)',...
    'IsoIntE (MJ/kg)','Critical Point','Triple Line')
title('Water')
scale = axis;
% axis([scale(1) scale(2) -1 25])
plotfixer
