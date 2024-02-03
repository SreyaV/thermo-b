
clf
fprintf('\nPlotting Exergy Paths...\n')

% Exergy Path
figure(1)
hold on
Xwater = mdot_w * [x0w x1w x2w x3w];
Xfuel  = [mdot_fuel * [x0f x1f x2f] mdot_mix*x3m];
Xair   = [mdot_air * [x0a x1a x2a] mdot_mix*x3m];
Xmix = mdot_mix * [x3m x4m x5m x6m];
plot([0,1,2,3],Xair/mdot_fuel/1e6,'-ob');
plot([0,1,2,3],Xfuel/mdot_fuel/1e6,'-or');
plot([3,4,5,6],Xmix/mdot_fuel/1e6,'-ok');
plot([0,1,2,3],Xwater/mdot_fuel/1e6,'-og');
xlabel('Station');
ylabel('Fuel-Specific Exergy MJ/kg-Fuel');
legend('Air','Fuel','Mix','Steam');
plotfixer
hold off

%Flow Exergy path
figure(2)
hold on
Xfwater = mdot_w * [xf0w xf1w xf2w xf3w];
Xffuel  = [mdot_fuel * [xf0f xf1f xf2f] mdot_mix*xf3m];
Xfair   = [mdot_air * [xf0a xf1a xf2a] mdot_mix*xf3m];
Xfmix = mdot_mix * [xf3m xf4m xf5m xf6m];
plot([0,1,2,3],Xfair/mdot_fuel/1e6,'-ob');
plot([0,1,2,3],Xffuel/mdot_fuel/1e6,'-or');
plot([3,4,5,6],Xfmix/mdot_fuel/1e6,'-ok');
plot([0,1,2,3],Xfwater/mdot_fuel/1e6,'-og');
xlabel('Station');
ylabel('Fuel-Specific Flow Exergy MJ/kg-Fuel');
legend('Air','Fuel','Mix','Steam');
plotfixer
hold off

% Exergy Balance Plot
loss_fuel_comp   = W_fuel_compressor/mdot_fuel + xf1f-xf2f;
loss_air_comp    = W_air_compressor/mdot_fuel + (xf1f-xf2a)*mdot_air/mdot_fuel;
loss_premixer    = xf2a*mdot_air/mdot_fuel + xf2f -xf3m*mdot_mix / mdot_fuel;
loss_combustor   = (xf3m-xf4m)*mdot_mix / mdot_fuel;
loss_gas_turbine = -W_gross_gas_turbine/mdot_fuel + (xf4m-xf5m)*mdot_mix / mdot_fuel;
loss_hrsg        = (xf5m - xf6m)*mdot_mix / mdot_fuel + (xf1w - xf2w)*mdot_w / mdot_fuel;
loss_feedp       = 0; 
loss_exhaust    = xf6m *mdot_mix / mdot_fuel;

exergy_losses = [loss_feedp, loss_fuel_comp, loss_air_comp, loss_premixer, loss_combustor, loss_gas_turbine, loss_hrsg, loss_exhaust];
figure(1)
bar(exergy_losses*100/sum(exergy_losses));
xlabel('Conversion Device');
ylabel('Exergy Loss(%)');
title(sprintf('LHV-Efficiency: %.2f',100*Eff_LHV),sprintf('Exergy-Efficiency: %.2f',100*Eff_exergy));
diff = 3;
text(6.5,45,'Device Key:');
text(6.5,42-0*diff,'1-Feed Pump');
text(6.5,42-1*diff,'2-Fuel Compressor');
text(6.5,42-2*diff,'3-Air Compressor');
text(6.5,42-3*diff,'4-Premixer');
text(6.5,42-4*diff,'5-Combustor');
text(6.5,42-5*diff,'6-GasTurbine');
text(6.5,42-6*diff,'7-HRSG');
text(6.5,42-7*diff,'8-Exhaust');

plotfixer
legend('off');