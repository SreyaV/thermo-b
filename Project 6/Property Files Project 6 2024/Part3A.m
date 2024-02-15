Tcrit_e = 33.19;               % K
Pcrit_e = 1.3152e6;            % Pa
rcrit_e = 14.936*2.0159;     % kg/m3

steps = 20;
delta = 200;

dT = Tcrit_e/delta;
dP = Pcrit_e/delta;
dr = rcrit_e/delta;

Tvals =  (Tcrit_e - steps*dT) : dT : (Tcrit_e + steps*dT);
Pvals = (Pcrit_e - steps*dP) : dP : (Pcrit_e + steps*dP);
rvals = (rcrit_e - steps*dr) : dr : (rcrit_e + steps*dr);

d1Pvals = zeros(length(Tvals), length(rvals));
d2Pvals = zeros(length(Tvals), length(rvals));

for i=1:length(Tvals)
    for j=1:length(rvals)
        d1Pvals(i, j) = dPdr_irT(ispecies,rvals(j),Tvals(i));
        d2Pvals(i, j) = d2Pdr2_irT(ispecies,rvals(j),Tvals(i));
    end
end

d1Pvals = d1Pvals ./ 1e3;
d2Pvals = d2Pvals ./ 1e3;

figure(1);
hold on;
[C1, h1] = contour(rvals, Tvals, d1Pvals, 'b');
clabel(C1, h1);
%colorbar;
title('First and Second Derivative of Pressure w.r.t. Density');
xlabel('Density (kg/m^3)');
ylabel('Temperature (K)');

% Plotting the second derivative d2Pvals
%figure(2);
[C2, h2] = contour(rvals, Tvals, d2Pvals, 'r');
clabel(C2, h2);
%colorbar;
%title('Second Derivative of Pressure w.r.t. Density (d^2P/dr^2)');
%xlabel('Density (kg/m^3)');
%ylabel('Temperature (K)');
plot(rcrit_e, Tcrit_e,'ko', 'MarkerSize', 10);
plot(comprcrit, compTcrit, 'p', 'MarkerSize', 15);


yline(Tcrit_e, '--')
xline(rcrit_e, '--')
legend('dP/d_{rho}', 'd2P/d_{rho}2', 'CP-Exp.', 'CP-Num.')
axis([27 33 30 36])

hold off;
plotfixer()
%% 

Tcrit_e = 33.19;               % K
rcrit_e = 14.936*2.0159;     % kg/m3

steps = 20;
delta = 1e-3;

dT = delta;
dr = delta;

Tvals =  32.7 : dT : 33.1;
rvals = 31.1 : dr : 31.5;
compTcrit = 0;
comprcrit = 0;

for i=1:length(Tvals)
    for j=1:length(rvals)
        d1Pval = dPdr_irT(ispecies,rvals(j),Tvals(i));
        d2Pval = d2Pdr2_irT(ispecies,rvals(j),Tvals(i));
        if abs(d1Pval)<0.01 & abs(d2Pval)<0.01
            compTcrit = Tvals(i);
            comprcrit = rvals(j);
            break
        end

    end
end

compTcrit
comprcrit
%
%% 
compTcrit = 32.93800;
comprcrit = 31.36000;
P_irT(ispecies, comprcrit, compTcrit)