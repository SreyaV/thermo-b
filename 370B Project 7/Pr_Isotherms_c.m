% Make a pressure-density diagram.  This file is useful for exploring how the P_crT
% function behaves in various regions (freezing line, P-r inflection point, etc.).
% This version is set up for the O2-N2-AR system and defaults to air.
% C.F. Edwards, 2/19/12 

addpath 'Fundamental Relation Files'
addpath 'Fundamental Relation Data'
addpath 'Mixture Models'
addpath 'Setup Files' 
addpath 'Property Files'
addpath 'Procedure Files'

clear all
format compact
fprintf('\n************************************************************\n')

% Set up the basic storage and load the FR files and mixture model.
% Set the number of components in mixture (N = 2 for binary, N = 3 for ternary).
N     = 3;
Setup_Air_Props;

% Set fixed-point values that are specific to air.  See Lemmon et al. 2000.
% The composition is: N2:O2:Ar = 0.7812:0.2096:0.0092
% Be sure to change these as needed if you use engineering air (79:21).
% And ignore these if you use an arbitrary composition.  They are for air.
Mair = 28.9586;                 % kg/kmol
Tmaxcondentherm = 132.6312;     % K
Pmaxcondentherm = 3.78502e6;    % Pa
rmaxcondentherm = 10.4477*Mair; % kg/m3
Tmaxcondenbar   = 132.6035;     % K
Pmaxcondenbar   = 3.7891e6;     % Pa
rmaxcondenbar   = 11.0948*Mair; % kg/m3
Tcritair        = 132.5306;     % K
Pcritair        = 3.7860e6;     % Pa
rcritair        = 11.8308*Mair; % kg/m3
% Bottom of dome conditions for air:
Tsolidair       = 59.75;        % K
Psolidair       = 5265;         % Pa
rsolidair       = 33.067*Mair;  % kg/m3
% Upper limit conditions for the model:
Tupper          = 870;          % K
Pupper          = 100e6;        % Pa
rupper          = rsolidair;
% Lower limit to stay just above solid air:
Tlower          = 60;           % K

% Set an air composition.
c = zeros(N,1);
if(N == 3)
    % Use ternary air.
    c(O2) = 0.2096;
    c(Ar) = 0.0092;
    c(N2) = 1 - c(O2) - c(Ar);
end
if(N == 2)
    % Use binary air.
    c(O2) = 0.21;
    c(N2) = 1 - c(O2);
end

% % Set an arbitrary mixture (or comment out).
% if(N == 3)
%     % Use a ternary mixture.
%     c(O2) = 0.25;
%     c(Ar) = 0.25;
%     c(N2) = 1 - c(O2) - c(Ar);
% end
% if(N == 2)
%     % Use a binary mixture.
%     c(O2) = 0.25;
%     c(N2) = 1 - c(O2);
% end

% % Get the inflection point for this composition.
[Tinfl rinfl] = Pr_Inflection_c(c)
Pinfl = P_crT(c,rinfl,Tinfl) 

% Start a P-r plot.  Get the spinodal lines first.
Tmax = Tinfl;
Tmin = Tlower;
steps = 100;
dT = (Tmax - Tmin)/steps;
i = 1;
Tplot   = zeros(steps+1);
rfsplot = zeros(steps+1);
rgsplot = zeros(steps+1);
Pfsplot = zeros(steps+1);
Pgsplot = zeros(steps+1);
for T=Tmin:dT:Tmax
    T
    % Start with the spinodal lines.
    rfs = Liquid_Spinodal_cT(c,T);
    Pfs = P_crT(c,rfs,T);
    rgs = Vapor_Spinodal_cT(c,T);    
    Pgs = P_crT(c,rgs,T);
        
    % Storage for plotting:
    Tplot(i)   = T;
    rfsplot(i) = rfs;
    rgsplot(i) = rgs;
    Pfsplot(i) = Pfs;
    Pgsplot(i) = Pgs;
    i = i+1;
end

% Put a few isotherms on the P-r plot.
Tlist = [...
    Tupper...
    0.4*(Tupper-Tinfl) + Tinfl ...
    0.1*(Tupper-Tinfl) + Tinfl ...
    2.00*(Tinfl-Tlower) + Tlower...
    1.50*(Tinfl-Tlower) + Tlower...
    1.25*(Tinfl-Tlower) + Tlower...
    1.05*(Tinfl-Tlower) + Tlower...
    1.00*(Tinfl-Tlower) + Tlower...
    0.95*(Tinfl-Tlower) + Tlower...
    0.9*(Tinfl-Tlower) + Tlower...
    0.75*(Tinfl-Tlower) + Tlower...
    0.5*(Tinfl-Tlower) + Tlower...
    0.25*(Tinfl-Tlower) + Tlower...
    0.1*(Tinfl-Tlower) + Tlower...
    Tlower...
    ]

rmax = rupper;
rmin = 1e-6;
steps = 1000;
dr = (rmax-rmin)/steps;
Pisotherm = zeros(length(Tlist),steps+1);
risotherm = zeros(length(Tlist),steps+1);
for j=1:1:length(Tlist)
    T = Tlist(j);
    i = 1;
    for r=rmin:dr:rmax
        Pisotherm(j,i) = P_crT(c,r,T);
        risotherm(j,i) = r;
        i = i+1;
    end
end

figure(1)
clf
hold on
plot(rinfl,Pinfl/1e6,'kd')
plot([0 rupper],[Psolidair Psolidair]/1e6,'k--')
for j=2:1:length(Tlist)-1
    plot(risotherm(j,:),Pisotherm(j,:)/1e6,'b')
end
plot(risotherm(1,:),Pisotherm(1,:)/1e6,'r')
plot(risotherm(length(Tlist),:),Pisotherm(length(Tlist),:)/1e6,'r')
plot(rgsplot,Pgsplot/1e6,'k.',rfsplot,Pfsplot/1e6,'k.')
plot(rinfl,Pinfl/1e6,'ko')
plot([rupper rupper],[-30 110],'r--')
plot([0 1000],[100 100],'r--')
legend('P-\rho Inflection','Solid Air','Location','SouthWest')
hold off
xlabel('Density (kg/m^3)')
ylabel('Pressure (MPa)')
for i=2:1:length(Tlist)-1
    text(750,1.9*(Pinfl/1e6)*(length(Tlist)-1.25*(i-1))/length(Tlist),...
        sprintf('       %.1f',Tlist(i)))
end
i = 1;
text(750,1.9*(Pinfl/1e6)*(length(Tlist)-1.25*(i-1))/length(Tlist),...
    sprintf('T = %.1f K',Tlist(i)),'Color','r')
i = length(Tlist);
text(750,1.9*(Pinfl/1e6)*(length(Tlist)-1.25*(i-1))/length(Tlist),...
    sprintf('       %.1f',Tlist(i)),'Color','r')
axis([0 rupper -0.5*Pinfl/1e6 2*Pinfl/1e6])
if(N == 3)
    title(sprintf('%.3f N2, %.3f O2, %.3f Ar',c(N2),c(O2),c(Ar)))
end
if(N == 2)
    title(sprintf('%.3f N2, %.3f O2',c(N2),c(O2)))
end
plotfixer
