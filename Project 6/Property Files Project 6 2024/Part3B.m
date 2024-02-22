Tcrit = 33.19;               % K
Pcrit = 1.3152e6;            % Pa

Ttrip  = 13.95;              % K
Ptrip  = 0.007199e6;         % Pa

steps = 500;
delta = 200;

dT = 0.01 ;%(Tcrit-Ttrip)/steps;
%dP = (Pcrit-Ptrip)/steps;


Tvals =  Ttrip+dT:dT:Tcrit-1.5;
%Pvals = Ptrip:dP:Pcrit;
Pvals = zeros(1, length(Tvals));

for i=1:length(Tvals)
    Tvals(i);
    temp = saturation_iT(ispecies, Tvals(i));
    Pvals(i) = temp(1);
end
%% 
Pvals = [Ptrip Pvals Pcrit];
Tvals= [Ttrip Tvals Tcrit];

%% 

nTcrit = 32.93800;
nrcrit = 31.36000;
nPcrit = 1.2838e6;

nTtrip = 13.950;
nPtrip = 7.671908e3;


%% 
clf
figure(1);
hold on;
xlabel('Temperature (K)')
ylabel('Pressure (MPa)')
axis([0 40 0 1.4])
plot(Tcrit, Pcrit/1e6, 'ko')
plot(nTcrit, nPcrit/1e6, 'rx')
plot(Ttrip, Ptrip/1e6, 'ko')
plot(nTtrip, nPtrip/1e6, 'rx')
plot(Tvals, Pvals./1e6)
legend('Experimental', 'Numerical')
text(15, 0.7, 'Liquid', 'HorizontalAlignment', 'center', 'FontSize', 30);
text(30, 0.2, 'Vapor', 'HorizontalAlignment', 'center', 'FontSize', 30);
title('P-T Plot for nH_2')
hold off;
plotfixer()

%% 
clf
figure(2);
hold on;
title("Dühring Plot for nH_2")
xlabel('-1/T (1/K)')
ylabel('ln(P) (/)')
%axis([13 35 1e-3 10])
plot(-1/Tcrit, log(Pcrit/1e6), 'ko')
plot(-1/nTcrit, log(nPcrit/1e6), 'rx')
plot(-1/Ttrip, log(Ptrip/1e6), 'ko')
plot(-1/nTtrip, log(nPtrip/1e6), 'rx')
inv_Tvals = zeros(1,length(Tvals));
for i=1:length(Tvals)
    inv_Tvals(i)=-1/Tvals(i);
end
plot(inv_Tvals, log(Pvals./1e6))
legend('Experimental', 'Numerical')
text(-0.06, -2, 'Liquid', 'HorizontalAlignment', 'center', 'FontSize', 30);
text(-0.04, -4, 'Vapor', 'HorizontalAlignment', 'center', 'FontSize', 30);
hold off;
plotfixer()

%% 

figure(3);
clf; % Clear figure to make sure we're starting fresh
semilogy(Tvals, Pvals./1e6, 'b-'); % Plot with logarithmic y-scale
hold on;

% Uncomment these if Tcrit and Pcrit, Ttrip and Ptrip are defined and you want to plot them
semilogy(Tcrit, Pcrit/1e6, 'ko'); % Plot critical point
semilogy(nTcrit, nPcrit/1e6, 'rx'); 
semilogy(Ttrip, Ptrip/1e6, 'ko'); % Plot triple point
semilogy(nTtrip, nPtrip/1e6, 'rx'); 
set(gca, 'XScale', 'log');
% Adjusting axis after plotting ensures the axis limits are correctly applied.
axis([13 35 1e-3 10]);

% Adding legend after all plots are made ensures it accurately reflects plotted data.
% Update legend entries according to the data you've plotted.
legend('','Experimental', 'Numerical')

title("Dühring Plot for nH_2");
xlabel('Temperature (K)');
ylabel('Pressure (MPa)');
text(20, 0.7, 'Liquid', 'HorizontalAlignment', 'center', 'FontSize', 30);
text(27, 0.2, 'Vapor', 'HorizontalAlignment', 'center', 'FontSize', 30);

% Hold off before applying plot adjustments like plotfixer to ensure they're applied correctly.
hold off;

% If plotfixer() adjusts plot properties, call it after setting up your plot completely.
plotfixer(); % Assuming plotfixer is a custom function to improve plot appearance


