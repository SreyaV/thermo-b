Tcrit = 33.19;               % K
Pcrit = 1.3152e6;            % Pa

Ttrip  = 13.95;              % K
Ptrip  = 0.007199e6;         % Pa

steps = 1000;
delta = 200;

dT = (Tcrit-Ttrip)/steps;
%dP = (Pcrit-Ptrip)/steps;

Tvals =  Ttrip+dT:dT:Tcrit-0.3;
%Pvals = Ptrip:dP:Pcrit;
Pvals = zeros(1, length(Tvals));

for i=1:length(Tvals)
    Tvals(i)
    temp = saturation_iT(ispecies, Tvals(i));
    Pvals(i) = temp(1);
end
%% 
Pvals = [Ptrip Pvals Pcrit];
Tvals= [Ttrip Tvals Tcrit];

%% 
clf
figure(1);
hold on;
xlabel('Temperature (K)')
ylabel('Pressure (MPa)')
axis([0 40 0 1.4])
plot(Tcrit, Pcrit/1e6, 'ko')
plot(Ttrip, Ptrip/1e6, 'ro')
plot(Tvals, Pvals./1e6)
legend('Critical Point', 'Triple Point')
hold off;
plotfixer()

%% 

figure(2);
hold on;
title("Dühring Plot")
xlabel('-1/T (1/K)')
ylabel('ln(P) (/)')
%axis([13 35 1e-3 10])
plot(-1/Tcrit, log(Pcrit/1e6), 'ko')
plot(-1/Ttrip, log(Ptrip/1e6), 'ro')
inv_Tvals = zeros(1,length(Tvals));
for i=1:length(Tvals)
    inv_Tvals(i)=-1/Tvals(i);
end
plot(inv_Tvals, log(Pvals./1e6))
legend('Critical Point', 'Triple Point')
hold off;
plotfixer()
