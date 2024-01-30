% Make T-s, u-s, h-s diagrams for water.
% C.F. Edwards, 10/30/12

clear all
format compact
fprintf('\n********************************************************\n')

environment = Solution('gri30.yaml','gri30');
N = nSpecies(environment);
iN2 = speciesIndex(environment, 'N2');
iO2 = speciesIndex(environment, 'O2');
iH2O = speciesIndex(environment, 'H2O');
iAr = speciesIndex(environment, 'AR');
iCO2 = speciesIndex(environment, 'CO2');
xenv = zeros(1, N);
xenv(iN2) = 0.757223;
xenv(iO2) = 0.202157;
xenv(iH2O) = 0.031208;
xenv(iAr) = 0.009015;
xenv(iCO2) = 0.000397;
set(environment, 'T', 300, 'P', oneatm, 'X', xenv);
T0 = 300;
P0 = oneatm;


% Set up a water object to work with in Cantera.
water = Solution('liquidvapor.yaml','water');

% Get the critical point props.
Tc = critTemperature(water)
Pc = critPressure(water)
set(water,'T',Tc,'P',Pc)
sc = entropy_mass(water);
uc = intEnergy_mass(water);
hc = enthalpy_mass(water);
xc = exergy_mass(water, 0, T0, P0, environment);

% Set the triple point props.  Use a value epsilon above Tt to get the
% bottom edge of the vapor dome.
Tt = 273.165
Pt = satPressure(water,Tt)
setState_Tsat(water,[Tt 0]);
sft = entropy_mass(water);
uft = intEnergy_mass(water);
hft = enthalpy_mass(water);
vft = 1/density(water);
xft = exergy_mass(water, 0, T0, P0, environment);
setState_Tsat(water,[Tt 1]);
sgt = entropy_mass(water);
ugt = intEnergy_mass(water);
hgt = enthalpy_mass(water);
vgt = 1/density(water);
xgt = exergy_mass(water, 0, T0, P0, environment);

% Set the limits for data curves.
Tmin = Tt;
Tmax = 1100;
Pmin = Pt;
Pmax = 10000*1e5;

% Make a vapor dome.
dT = (Tc-Tmin)/150;
i = 1;
%setState_Tsat(water,[Tmin 0]);
for T=Tmin:dT:Tc-dT    % Stop short of the critical point.
    Tsatline(i) = T;
    setState_Tsat(water,[T 0]);
    sliqline(i) = entropy_mass(water);
    uliqline(i) = intEnergy_mass(water);
    hliqline(i) = enthalpy_mass(water);
    xliqline(i) = exergy_mass(water, 0, T0, P0, environment);
    setState_Tsat(water,[T 1]);
    svapline(i) = entropy_mass(water);
    uvapline(i) = intEnergy_mass(water);
    hvapline(i) = enthalpy_mass(water);
    xvapline(i) = exergy_mass(water, 0, T0, P0, environment);
    i = i+1;
end
Tsatline(i) = Tc;   % Add the critical point now.
sliqline(i) = sc;
uliqline(i) = uc;
hliqline(i) = hc;
xliqline(i) = xc;
svapline(i) = sc;
uvapline(i) = uc;
hvapline(i) = hc;
xvapline(i) = xc;

% Start a set of isobaric curves.
Plist = [0.01 0.1 1 10 100];  % Pressure in bar
for j=1:1:length(Plist)
    water = Solution('liquidvapor.yaml','water');
    P = Plist(j)*100000

    % Do the compressed liquid side.
    setState_Psat(water,[P 0]);
    Tsat = temperature(water);
    dT = (Tsat - Tmin)/100;
    i = 1;
    for T=Tmin:dT:Tsat-dT  % Stop short of saturation.
        Tpresline(j,i) = T;
        set(water,'T',T,'P',P);
        spresline(j,i) = entropy_mass(water);
        upresline(j,i) = intEnergy_mass(water);
        hpresline(j,i) = enthalpy_mass(water);
        xpresline(j,i) = exergy_mass(water, 0, T0, P0, environment);
        i = i+1;
    end
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    setState_Psat(water,[P 0]);
    spresline(j,i) = entropy_mass(water);
    upresline(j,i) = intEnergy_mass(water);
    hpresline(j,i) = enthalpy_mass(water);
    xpresline(j, i) = exergy_mass(water, 0, T0, P0, environment);
    i = i+1;

    % Now go across the dome.
    dq = 1/100;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Tpresline(j,i) = Tsat;
        setState_Psat(water,[P q]);
        spresline(j,i) = entropy_mass(water);
        upresline(j,i) = intEnergy_mass(water);
        hpresline(j,i) = enthalpy_mass(water);
        xpresline(j, i) = exergy_mass(water, 0, T0, P0, environment);
        i = i+1;
    end
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    setState_Psat(water,[P 1]);
    spresline(j,i) = entropy_mass(water);
    upresline(j,i) = intEnergy_mass(water);
    hpresline(j,i) = enthalpy_mass(water);
    xpresline(j, i) = exergy_mass(water, 0, T0, P0, environment);
    i = i+1;

    % Do the vapor side.
    dT = (Tmax - Tsat)/100;
    for T=Tsat+dT:dT:Tmax  % Start just above saturation.
        try
            set(water,'T',T,'P',P);
            Tpresline(j,i) = T;
            spresline(j,i) = entropy_mass(water);
            upresline(j,i) = intEnergy_mass(water);
            hpresline(j,i) = enthalpy_mass(water);
            xpresline(j, i) = exergy_mass(water, 0, T0, P0, environment);
        catch
            disp('Trouble finding isobar...')
            water = Solution('liquidvapor.xml','water');
            Tpresline(j,i) = Tpresline(j,i-1);
            spresline(j,i) = spresline(j,i-1);
            upresline(j,i) = upresline(j,i-1);
            hpresline(j,i) = hpresline(j,i-1);
            xpresline(j, i) = exergy_mass(water, 0, T0, P0, environment);
        end
        i = i+1;
    end
end

% Add an isobar above the critical pressure.
P = 1000;           % In bar
Plist = [Plist P];  % Add to the list
P = P*100000;       % In Pascals
j = length(Plist);
dT = (Tmax - Tmin)/300;
i = 1;
for T=Tmin:dT:Tmax  % Stop short of saturation.
    Tpresline(j,i) = T;
    set(water,'T',T,'P',P);
    spresline(j,i) = entropy_mass(water);
    upresline(j,i) = intEnergy_mass(water);
    hpresline(j,i) = enthalpy_mass(water);
    xpresline(j, i) = exergy_mass(water, 0, T0, P0, environment);
    i = i+1;
end

% Add an isobar above the critical pressure.
P = 10000;           % In bar
Plist = [Plist P];  % Add to the list
P = P*100000;       % In Pascals
j = length(Plist);
dT = (Tmax - Tmin)/300;
i = 1;
for T=Tmax:-dT:Tmin  % Stop short of saturation.
    set(water,'T',T,'P',P);
    if T > 340
        Tpresline(j,i) = T;
        spresline(j,i) = entropy_mass(water);
        upresline(j,i) = intEnergy_mass(water);
        hpresline(j,i) = enthalpy_mass(water);
        xpresline(j, i) = exergy_mass(water, 0, T0, P0, environment);
    else
        % Can't do liquids when too dense.
        Tpresline(j,i) = Tpresline(j,i-1);
        spresline(j,i) = spresline(j,i-1);
        upresline(j,i) = upresline(j,i-1);
        hpresline(j,i) = hpresline(j,i-1);
        xpresline(j, i) = xpresline(j,i-1);
    end
    i = i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isotherms.
Tlist = [274 400 500 600 700 800 900 1000];
for j=1:1:3    % Do the part of list below the critical point. 
    water = Solution('liquidvapor.yaml','water');
    T = Tlist(j)

    % Do the compressed liquid side.
    setState_Tsat(water,[T 0]);
    Psat = pressure(water);
    logdP = (log(Pmax) - log(Psat))/100;
    i = 1;
    for logP=log(Pmax):-logdP:log(Psat)+logdP  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(water,'T',T,'P',P);
        stempline(j,i) = entropy_mass(water);
        utempline(j,i) = intEnergy_mass(water);
        htempline(j,i) = enthalpy_mass(water);
        xtempline(j,i) = exergy_mass(water, 0, T0, P0, environment);
        i = i+1;
    end
    Ttempline(j,i) = T;   % Add the saturation point now.
    % set(water,'T',T,'P',Psat);
    setState_Tsat(water,[T 0]);
    stempline(j,i) = entropy_mass(water);
    utempline(j,i) = intEnergy_mass(water);
    htempline(j,i) = enthalpy_mass(water);
    xtempline(j,i) = exergy_mass(water, 0, T0, P0, environment);
    i = i+1;

    % Now go across the dome.
    dq = 1/100;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Ttempline(j,i) = T;
        setState_Psat(water,[Psat q]);
        stempline(j,i) = entropy_mass(water);
        utempline(j,i) = intEnergy_mass(water);
        htempline(j,i) = enthalpy_mass(water);
        xtempline(j,i) = exergy_mass(water, 0, T0, P0, environment);
        i = i+1;
    end
    Ttempline(j,i) = T; % Add the saturation point now.
    setState_Psat(water,[Psat 1]);
    stempline(j,i) = entropy_mass(water);
    utempline(j,i) = intEnergy_mass(water);
    htempline(j,i) = enthalpy_mass(water);
    xtempline(j,i) = exergy_mass(water, 0, T0, P0, environment);
    i = i+1;

    % Do the vapor side.
    logdP = (log(Psat) - log(Pmin))/100;
    for logP=log(Psat)-logdP:-logdP:log(Pmin)  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(water,'T',T,'P',P);
        stempline(j,i) = entropy_mass(water);
        utempline(j,i) = intEnergy_mass(water);
        htempline(j,i) = enthalpy_mass(water);
        xtempline(j,i) = exergy_mass(water, 0, T0, P0, environment);
        i = i+1;
    end
end

% Add isotherms above the critical temperature.
for j=j+1:1:length(Tlist)
    water = Solution('liquidvapor.yaml','water');
    T = Tlist(j)

    logdP = (log(Pmax) - log(Pmin))/300;
    i = 1;
    for logP=log(Pmax):-logdP:log(Pmin)  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(water,'T',T,'P',P);
        stempline(j,i) = entropy_mass(water);
        utempline(j,i) = intEnergy_mass(water);
        htempline(j,i) = enthalpy_mass(water);
        xtempline(j,i) = exergy_mass(water, 0, T0, P0, environment);
        i = i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isochores.
set(water,'T',Tmin+.01,'P',Pmin);
vmin = 1/density(water);
vmax = vgt;
set(water,'T',Tmax,'V',vmin);
vlist = [0.001 0.01 0.1 1 10 100];
for j=1:1:length(vlist)
    v = vlist(j)
    dT = (Tmax - Tmin)/300;
    i = 1;
    for T=Tmin:dT:Tmax  % Stop short of saturation.
        Tvoluline(j,i) = T;
        set(water,'T',T,'V',v);
        svoluline(j,i) = entropy_mass(water);
        uvoluline(j,i) = intEnergy_mass(water);
        hvoluline(j,i) = enthalpy_mass(water);
        xvoluline(j,i) = exergy_mass(water, 0, T0, P0, environment);
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
            xerguline(j,i) = exergy_mass(water, 0, T0, P0, environment);
            if(Terguline(j,i) > Tmax)
                Terguline(j,i) = Terguline(j,i-1);
                serguline(j,i) = serguline(j,i-1);
                xerguline(j,i) = exergy_mass(water, 0, T0, P0, environment);
            end
        catch
            disp('Trouble finding isoerg...')
            water = Solution('liquidvapor.yaml','water');
            Terguline(j,i) = Terguline(j,i-1);
            serguline(j,i) = serguline(j,i-1);
            xerguline(j,i) = exergy_mass(water, 0, T0, P0, environment);
        end
        i = i+1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isenthalps.
hlist = [-13 -13.25 -13.5 -14 -14.5 -15 -15.5]*1e6;
Pmax  = 10000e5;
Pmin = .01e5;
for j=1:1:length(hlist)
    h = hlist(j)
    lPhigh = log(Pmax);
    lPlow  = log(Pmin);
    ldP = (lPhigh - lPlow)/300;
    i = 1;
    for lP=lPlow:ldP:lPhigh  % Stop short of saturation.
        P = exp(lP);
        try
            setState_HP(water,[h,P]);
            Tenthline(j,i) = temperature(water);
            senthline(j,i) = entropy_mass(water);
            xenthline(j,i) = exergy_mass(water, 0, T0, P0, environment);
            if(Tenthline(j,i) < Tmin)
                Tenthline(j,i) = Tenthline(j,i-1);
                senthline(j,i) = senthline(j,i-1);
                xenthline(j,i) = xenthline(j,i-1);
            end
        catch
            disp('Trouble finding isenthalp...')
            water = Solution('liquidvapor.yaml','water');
            Tenthline(j,i) = Tenthline(j,i-1);
            senthline(j,i) = senthline(j,i-1);
            xenthline(j,i) = xenthline(j,i-1);
        end
        i = i+1;
    end
end

% xpresline = xpresline ./1e6;
% xvoluline = xvoluline ./1e6
% xerguline = xerguline ./ 1e6;
% xenthline = xenthline ./ 1e6;
% xft = xft/1e6;
% xgt = xgt /1e6;


%% 

% Make the plots.
figure(1)
clf
hold on
i = 1;
plot(spresline(i,:)/1000,xpresline(i,:)/1e6,'b')
plot(svoluline(i,:)/1000,xvoluline(i,:)/1e6,'-','Color',[0 .6 0])
%plot(serguline(i,:)/1000,xerguline(i,:)/1e6,'m')
plot(stempline(i,:)/1000,xtempline(i,:)/1e6,'r')
plot(sc/1000,xc,'kd')
plot([sft/1000 sgt/1000],[xft/1e6 xgt/1e6],'ko--')
for i=1:1:length(Plist)
    plot(spresline(i,:)/1000,xpresline(i,:)/1e6,'b')
end
for i=1:1:length(vlist)
    plot(svoluline(i,:)/1000,xvoluline(i,:)/1e6,'-','Color',[0 .6 0])
end
for i=1:1:length(ulist)
    if i~=3
        plot(serguline(i,10:290)/1000,xerguline(i,10:290)/1e6,'m')
    end
end
for i=1:1:length(Tlist)
    plot(stempline(i,:)/1000,xtempline(i,:)/1e6,'r')
end
plot(sliqline/1000,xliqline/1e6,'k')
plot(svapline/1000,xvapline/1e6,'k')

plot([sft/1000 sgt/1000],[xft/1e6 xgt/1e6],'ko--')
hold off
xlabel('Specific Entropy (kJ/kg-K)')
ylabel('Specific Exergy (MJ/kg)')
for i=1:1:length(Plist)
    text(spresline(i,150)/1000,xpresline(i,150)/1e6,num2str(Plist(i)),'Color','b')
end
for i=1:1:length(vlist)
    text(svoluline(i,150)/1000,xvoluline(i,150)/1e6,num2str(vlist(i)),'Color',[0 .6 0])
end
for i=1:1:length(ulist)
    text(serguline(i,1)/1000,xerguline(i,1)/1e6,num2str(ulist(i)/1e6),'Color','m')
end

for i=1:1:length(Tlist)
    if Tlist(i) < 600
        text(stempline(i,150)/1000+3.6+i*.55,xtempline(i,150)/1e6+48*i/6,num2str(Tlist(i)),'Color','r')

    else
        text(stempline(i,150)/1000+2.7+i/15,xtempline(i,150)/1e6+78*i/8,num2str(Tlist(i)),'Color','r')
    end
end

% for i=1:1:length(Tlist)
%     text(stempline(i,1)/1000,xtempline(i,1)/1e6,num2str(Tlist(i)),'Color','r')
% end

legend('Isobar (bar)','Isochore (m^3/kg)','IsoIntE (MJ/kg)','Isotherm (K)','Critical Point','Triple Line')
title('Water')
scale = axis;
axis([4 14 0 1.6])
plotfixer_v2()


% figure(2)
% clf
% hold on
% i = 1;
% plot(spresline(i,:)/1000,upresline(i,:)/1e6,'b')
% plot(stempline(i,:)/1000,utempline(i,:)/1e6,'r')
% plot(svoluline(i,:)/1000,uvoluline(i,:)/1e6,'g')
% plot(sc/1000,uc/1e6,'kd')
% plot([sft/1000 sgt/1000],[uft/1e6 ugt/1e6],'ko--')
% for i=1:1:length(Plist)
%     plot(spresline(i,:)/1000,upresline(i,:)/1e6,'b')
% end
% for i=1:1:length(Tlist)
%     plot(stempline(i,:)/1000,utempline(i,:)/1e6,'r')
% end
% for i=1:1:length(vlist)
%     plot(svoluline(i,:)/1000,uvoluline(i,:)/1e6,'g')
% end
% plot(sliqline/1000,uliqline/1e6,'k')
% plot(svapline/1000,uvapline/1e6,'k')
% plot([sft/1000 sgt/1000],[uft/1e6 ugt/1e6],'ko--')
% hold off
% xlabel('Specific Entropy (kJ/kg-K)')
% ylabel('Specific Internal Energy (MJ/kg)')
% for i=1:1:length(Plist)
%     text(spresline(i,150)/1000,upresline(i,150)/1e6,...
%         num2str(Plist(i)),'Color','b')
% end
% for i=1:1:length(Tlist)
%     text(stempline(i,150)/1000,utempline(i,150)/1e6,...
%         num2str(Tlist(i)),'Color','r')
% end
% for i=1:1:length(vlist)
%     text(svoluline(i,150)/1000,uvoluline(i,150)/1e6,...
%         num2str(vlist(i)),'Color',[0 .6 0])
% end
% legend('Isobar (bar)','Isotherm (K)','Isochore (m^3/kg)',...
%     'Critical Point','Triple Line')
% title('Water')
% plotfixer
% 
% 
% figure(3)
% clf
% hold on
% i = 1;
% plot(spresline(i,:)/1000,hpresline(i,:)/1e6,'b')
% plot(stempline(i,:)/1000,htempline(i,:)/1e6,'r')
% plot(svoluline(i,:)/1000,hvoluline(i,:)/1e6,'g')
% plot(sc/1000,hc/1e6,'kd')
% plot([sft/1000 sgt/1000],[hft/1e6 hgt/1e6],'ko--')
% for i=1:1:length(Plist)
%     plot(spresline(i,:)/1000,hpresline(i,:)/1e6,'b')
% end
% for i=1:1:length(Tlist)
%     plot(stempline(i,:)/1000,htempline(i,:)/1e6,'r')
% end
% for i=1:1:length(vlist)
%     plot(svoluline(i,:)/1000,hvoluline(i,:)/1e6,'g')
% end
% plot(sliqline/1000,hliqline/1e6,'k')
% plot(svapline/1000,hvapline/1e6,'k')
% plot([sft/1000 sgt/1000],[hft/1e6 hgt/1e6],'ko--')
% hold off
% scale = axis;
% axis([scale(1) scale(2) scale(3) -11.5])
% xlabel('Specific Entropy (kJ/kg-K)')
% ylabel('Specific Enthalpy (MJ/kg)')
% for i=1:1:length(Plist)
%     text(spresline(i,300)/1000,hpresline(i,300)/1e6,...
%         num2str(Plist(i)),'Color','b')
% end
% for i=1:1:length(Tlist)
%     text(stempline(i,150)/1000,htempline(i,150)/1e6,...
%         num2str(Tlist(i)),'Color','r')
% end
% for i=1:1:length(vlist)
%     text(svoluline(i,200)/1000,hvoluline(i,200)/1e6,...
%         num2str(vlist(i)),'Color',[0 .6 0])
% end
% legend('Isobar (bar)','Isotherm (K)','Isochore (m^3/kg)',...
%     'Critical Point','Triple Line')
% title('Water')
% plotfixer

