%% Plots
% input feed composition required
hold on
plot(quality_data, x_data(:, N2), 'o-g', 'MarkerSize', 12)
plot(quality_data, x_data(:, O2), 'o-b', 'MarkerSize', 12)
plot(quality_data, x_data(:, Ar)*10, 'o-m', 'MarkerSize', 12)
plot([0 1], [0 0], '-k')
plot([0 1], [0 0], '--k')
plot([0 1], [0 0], '+--k')
plot(quality_data, y_data(:, N2), '+--g')
plot(quality_data, y_data(:, O2), '+--b')
plot(quality_data, y_data(:, Ar)*10, '+--m')
plot([0 1], [x_out(N2) x_out(N2)], '--g')
plot([0 1], [x_out(O2) x_out(O2)], '--b')
plot([0 1], [x_out(Ar) x_out(Ar)]*10, '--m')
legend('Nitrogen', 'Oxygen', 'Argon*10', 'Inlet', 'Product', 'Reflux')
ylabel('Mole Fraction')
xlabel('Reboiler Outlet Quality (mass)')
plotfixer()
%%
hold on
plot(quality_data, T_data(:, 1), 'o-b', 'MarkerSize', 12)
plot(quality_data, T_data(:, 2), '+--k', 'MarkerSize', 12)
ylabel('Temperature(K)')
xlabel('Reboiler Outlet Quality (mass)')
legend('Inlet', 'Outlet')
axis([0 1 87 89.5])
plotfixer()
%%
hold on
plot(quality_data, Q_data/1000, 'o-r', 'MarkerSize', 12)
ylabel('Feed-Specific Heat Rate (kJ/kg)')
xlabel('Reboiler Outlet Quality (mass)')
plotfixer()