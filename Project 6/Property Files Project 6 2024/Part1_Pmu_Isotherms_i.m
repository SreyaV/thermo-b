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
    ]

% Set limits and step size.
rmax = rupper_i(ispecies);
rmin = rgtrip_i(ispecies)/2;
steps = 2000;
dr = (rmax-rmin)/steps;

% Preallocate storage...
Prisotherm  = zeros(length(Tlist),steps+1);
risotherm   = zeros(length(Tlist),steps+1);
murisotherm = zeros(length(Tlist),steps+1);

for j=1:length(Tlist)
    T = Tlist(j)
    i = 1;
    
    % Spinodals computation
    left_spin  = spinodalL(ispecies,T);
    right_spin = spinodalV(ispecies,T);
 
    for r=rmin:dr:rmax+dr
        
        % Outside spinodals
        if r>right_spin || r<left_spin
            Prisotherm(j,i) = P_irT(ispecies,r,T);
            murisotherm(j, i) = mu_irT(ispecies, r, T);
            risotherm(j,i) = r;
            i = i+1;
  
        else
            Prisotherm(j,i) =NaN;
            murisotherm(j, i) = NaN;
            risotherm(j,i) = NaN;
            i = i+1;
        end

    end
end

%% Saturation line

disp('Generating saturation lines for Temperatures (K)...')
steps = 50;
dT = (-Ttrip_i(ispecies)+Tcrit_i(ispecies))/steps;
Tlist_sat = Ttrip_i(ispecies):dT:Tcrit_i(ispecies);

for j=1:length(Tlist_sat)
    
    T = Tlist_sat(j)
    [Psat rf rg] = Saturation_iT(ispecies,T);
    P_sat(j) = Psat;
    Mu_sat(j) = mu_irT(ispecies,rf,T);

end

%%
figure(1)
clf

% Put the isotherms on.
for j=1:1:length(Tlist)
    plot(Prisotherm(j,:)/1e6, murisotherm(j,:)/1e6,'b')
    hold on
end

% Add saturation line
plot(P_sat*1e-6,Mu_sat/1e6,'r')

hold off

ylabel('Chemical Potential (MJ/kmol)')
xlabel('Pressure (MPa)')

% Add some simple temperature labels.
for i=1:1:length(Tlist)
    text(1.7,-0.7/length(Tlist)*i,...
        num2str(Tlist(i)))
end
text(1.7,0.0,'T (K) =')

axis([-0.2 2 -0.8 0.2])
% Gussy up the plot a little.
plotfixer