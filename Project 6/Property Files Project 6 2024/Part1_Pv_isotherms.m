% Make a pressure-density isotherm diagram.  
% This file is useful for exploring how the P_rT function behaves in various 
% regions (triple line, critical point, etc.).
% C.F. Edwards, 2-11-12 

% Provide access to support files via the Matlab path.
addpath 'Fundamental_Relation_Files' 
addpath 'Fundamental_Relation_Data'
addpath 'Setup_Files' 
addpath 'Property_Files' 

% Clean up and get ready to go.
clear all
format compact
fprintf('\n************************************************************\n')

% Set up the basic storage and load the FR files.
Setup_Props_i;

% Set which of the loaded species you want to work with.  You might want to
% change this around to see what the plots look like for other species.
ispecies = nH2

% Put a few isotherms on the P-r plot.
Tlist = [...
    Ttrip_i(ispecies)...
    0.1*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.25*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.5*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.75*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.9*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.95*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.00*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.05*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.25*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    1.50*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    2.0*(Tcrit_i(ispecies)-Ttrip_i(ispecies)) + Ttrip_i(ispecies)...
    0.1*(Tupper_i(ispecies)-Tcrit_i(ispecies)) + Tcrit_i(ispecies) ...
    0.4*(Tupper_i(ispecies)-Tcrit_i(ispecies)) + Tcrit_i(ispecies) ...
    Tupper_i(ispecies)...
    ]

% Set limits and step size.
rmax = rupper_i(ispecies);
rmin = rgtrip_i(ispecies)/2;
steps = 2000;
dr = (rmax-rmin)/steps;

% Preallocate storage...
Prisotherm  = zeros(length(Tlist),steps+1);
risotherm   = zeros(length(Tlist),steps+1);
Vrisotherm = zeros(length(Tlist),steps+1);

for j=1:1:length(Tlist)

    T = Tlist(j)
    i = 1;

    % Spinodals computation
    left_spin  = spinodalL(ispecies,T);
    right_spin = spinodalV(ispecies,T);

    for r=rmin:dr:rmax+dr

        % Outside spinodals
        if r>right_spin || r<left_spin
            Prisotherm(j,i) = P_irT(ispecies,r,T);
            risotherm(j,i) = r;
            Vrisotherm(j,i) = 1/r;
            i = i+1;

        else
            Prisotherm(j,i) =NaN;
            risotherm(j,i) = NaN;
            Vrisotherm(j,i) = NaN;
            i = i+1;
        end
    end
end

 %% Saturation line

disp('Generating saturation lines for Temperatures (K)...')
steps = 50;
dT = (-Ttrip_i(ispecies)+Tcrit_i(ispecies)*0.99)/steps;
Tlist_sat = Ttrip_i(ispecies):dT:Tcrit_i(ispecies)*0.99;
N=length(Tlist_sat);

for j=1:length(Tlist_sat)
    
    
    T = Tlist_sat(j)
    try
        [Psat rf rg] = Saturation_iT(ispecies,T);
        P_sat_f(j) = Psat;
        P_sat_g(N+1-j)  = Psat;
        V_sat_f(j) = 1/rf;
        V_sat_g(N+1-j) = 1/rg;
    catch
        P_sat_f(j) = nan;
        P_sat_g(N+1-j)  = nan;
        V_sat_f(j) = nan;
        V_sat_g(N+1-j) = nan;
    end
end


%% Triple Line

Vtrip = 1/rmax :0.01:1/rmin;
for i=1:length(Vtrip)
    Ptrip(i)=Ptrip_i(ispecies);
end

%%
figure(1)
clf

semilogx(rcrit_i(ispecies),Pcrit_i(ispecies),'kd')
hold on

for j=1:1:length(Tlist)
    plot(Vrisotherm(j,:),Prisotherm(j,:)/1e6,'b')
end

% Add saturation line
Pcrit_num = 1.2838e6; %[Pa]
P_sat = [P_sat_f Pcrit_num P_sat_g];
V_sat = [V_sat_f 1/rcrit_i(ispecies) V_sat_g];
plot(V_sat,P_sat/1e6,'r')

% Add triple line
plot(Vtrip, Ptrip/1e6,'r')

hold off
xlabel('Specific Volume (kg/m^3)')
ylabel('Pressure (MPa)')


% Add some simple temperature labels.
for i=1:1:length(Tlist)
    text(5,1.7/length(Tlist)*i,...
        num2str(Tlist(i)))
end
text(4,1.8,'T (K) =')

axis([0.01 10 -0.2 2])

% Gussy up the plot a little.
plotfixer
