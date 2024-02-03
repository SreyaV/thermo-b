
clf
fprintf('\nPlotting Exergy Paths...\n')

% % Exergy Path
% figure(1)
% hold on
% Xwater = mdot_w * [xw1 xw2 xw3 xw4 xw5 xw6 xw7 xw8];
% Xfuel  = [mdot_fuel * [xf1 xf2] mdot_mix*xm3];
% Xair   = [mdot_air * [xa1 xa2] mdot_mix*xm3];
% Xmix = mdot_mix * [xm3 xm4 xm5 xm6 xm7 xm8];
% plot([1,2,3],Xair/mdot_fuel/1e6,'-ob');
% plot([1,2,3],Xfuel/mdot_fuel/1e6,'-or');
% plot([3,4,5,6,7,8],Xmix/mdot_fuel/1e6,'-ok');
% plot([1,2,3,4,5,6,7,8],Xwater/mdot_fuel/1e6,'-og');
% xlabel('Station');
% ylabel('Fuel-Specific Exergy MJ/kg-Fuel');
% legend('Air','Fuel','Mix','Steam');
% plotfixer
% hold off
% 
% %Flow Exergy path
% figure(2)
% hold on
% Xfwater = mdot_w * [xfw1 xfw2 xfw3 xfw4 xfw5 xfw6 xfw7 xfw8];
% Xffuel  = [mdot_fuel * [xff1 xff2] mdot_mix*xfm3];
% Xfair   = [mdot_air * [xfa1 xfa2] mdot_mix*xfm3];
% Xfmix = mdot_mix * [xfm3 xfm4 xfm5 xfm6 xfm7 xfm8];
% plot([1,2,3],Xfair/mdot_fuel/1e6,'-ob');
% plot([1,2,3],Xffuel/mdot_fuel/1e6,'-or');
% plot([3,4,5,6,7,8],Xfmix/mdot_fuel/1e6,'-ok');
% plot([1,2,3,4,5,6,7,8],Xfwater/mdot_fuel/1e6,'-og');
% xlabel('Station');
% ylabel('Fuel-Specific Flow Exergy MJ/kg-Fuel');
% legend('Air','Fuel','Mix','Steam');
% plotfixer
% hold off

% Exergy Balance Plot
loss_fuel_comp = W_fuel_compressor/mdot_fuel + xff1-xff2;
loss_air_comp = W_air_compressor/mdot_fuel + (xfa1-xfa2)*mdot_air/mdot_fuel;
loss_premixer = xfa2*mdot_air/mdot_fuel + xff2 -xfm3*mdot_mix / mdot_fuel;
loss_combustor = (xfm3-xfm4)*mdot_mix / mdot_fuel;
loss_gas_turbine= -W_gross_gas_turbine/mdot_fuel + (xfm4-xfm5)*mdot_mix / mdot_fuel;
loss_hrsg = (xfm5 - xfm8)*mdot_mix / mdot_fuel + (xfw2 - xfw5)*mdot_w / mdot_fuel;
loss_feedp = 0; %W_feed_pump/mdot_fuel + (xfw1 - xfw2)*mdot_w / mdot_fuel;
loss_steam_turbine = -W_gross_steam_turbine/mdot_fuel +(xfw5 - xfw6)*mdot_w / mdot_fuel;
loss_condenser= (xfw6 - xfw7)*mdot_w / mdot_fuel;
loss_cond_pump= W_cond_pump/mdot_fuel + (xfw7 - xfw8)*mdot_w / mdot_fuel;
loss_cold_storage = 0;
loss_flue_gas = xfm8 *mdot_mix / mdot_fuel;
exergy_losses = [loss_fuel_comp, loss_air_comp, loss_premixer, loss_combustor, loss_gas_turbine, loss_hrsg, loss_feedp, loss_steam_turbine, loss_condenser, loss_cond_pump, loss_cold_storage, loss_flue_gas];
figure(1)
bar(exergy_losses*100/sum(exergy_losses));
xlabel('Conversion Device');
ylabel('Eexergy Loss(%)');
title(sprintf('LHV-Efficiency: %.2f',100*Eff_LHV),sprintf('Exergy-Efficiency: %.2f',100*Eff_exergy));
diff = 3;
text(6.5,55,'Device Key:');
text(6.5,55-1*diff,'1-Fuel Compressor');
text(6.5,55-2*diff,'2-Air Compressor');
text(6.5,55-3*diff,'3-Premixer');
text(6.5,55-4*diff,'4-Combustor');
text(6.5,55-5*diff,'5-GasTurbine');
text(6.5,55-6*diff,'6-HRSG');
text(6.5,55-7*diff,'7-Feed Pump');
text(6.5,55-8*diff,'8-Steam Turbine');
text(6.5,55-9*diff,'9-Condenser');
text(6.5,55-10*diff,'10-Condenser Pump');
text(6.5,55-11*diff,'11-Cold Storage');
text(6.5,55-12*diff,'12-Flue gas');
plotfixer
legend('off');