hold on
Tc = critTemperature(water);
Pc = critPressure(water);
set(water,'T',Tc,'P',Pc);
hc = enthalpy_mass(water);
dT = (Tc-280-.1)/100;
i = 1;
for T=280:dT:Tc
    Tsatline(i) = T;
    setState_Tsat(water,[T 0]);
    hliqline(i) = enthalpy_mass(water);
    setState_Tsat(water,[T 1]);
    hvapline(i) = enthalpy_mass(water);
    i = i+1;
end
% Add the critical point at the top.
Tsatline(i) = Tc;
hliqline(i) = hc;
hvapline(i) = hc;

Hoffset_water = hw5*mdot_w;
Hoffset_mix = hm5*mdot_mix;

% Make T-h plot.

Twater = [Tw1 Tw2 Tw3 Tw4 Tw5 Tw6 Tw7 Tw8];
Hwater = mdot_w * [hw1 hw2 hw3 hw4 hw5 hw6 hw7 hw8];
Tfuel  = [Tf1 Tf2];
Hfuel  = mdot_fuel * [hf1 hf2];
Tair   = [Ta1 Ta2];
Hair   = mdot_air * [ha1 ha2];
Tmix = [Tm3 Tm4 Tm5 Tm6 Tm7 Tm8];
Hmix = mdot_mix * [hm3 hm4 hm5 hm6 hm7 hm8];

plot((Hair-Hoffset_mix)/mdot_fuel/1e6,Tair,'bo--')
plot((Hfuel-Hoffset_mix)/mdot_fuel/1e6,Tfuel,'ro--')
plot((Hmix-Hoffset_mix)/mdot_fuel/1e6,Tmix,'ko--')
plot((Hwater-Hoffset_water)/mdot_fuel/1e6,Twater,'ko--')


plot(([hf2*mdot_fuel hm3*mdot_mix]-Hoffset_mix)/mdot_fuel/1e6,[Tf2 Tm3],'--')
plot(([ha2*mdot_air hm3*mdot_mix]-Hoffset_mix)/mdot_fuel/1e6,[Ta2 Tm3],'--')

plot((mdot_w*hliqline - Hoffset_water)/mdot_fuel/1e6,Tsatline,'Color',[0 .5 0])
plot((mdot_w*hvapline - Hoffset_water)/mdot_fuel/1e6,Tsatline,'Color',[0 .5 0])

text((Hair(1)-Hoffset_mix)/mdot_fuel/1e6,Tair(1),'a,1')
text((Hfuel(1)-Hoffset_mix)/mdot_fuel/1e6,Tfuel(1),'f,1')
text((Hair(2)-Hoffset_mix)/mdot_fuel/1e6,Tair(2),'a,2')
text((Hfuel(2)-Hoffset_mix)/mdot_fuel/1e6,Tfuel(2),'f,2')
text((Hmix(1)-Hoffset_mix)/mdot_fuel/1e6,Tmix(1),'m,3')
text((Hmix(2)-Hoffset_mix)/mdot_fuel/1e6,Tmix(2),'m,4')
text((Hmix(3)-Hoffset_mix)/mdot_fuel/1e6,Tmix(3),'m,5')
text((Hmix(4)-Hoffset_mix)/mdot_fuel/1e6,Tmix(4),'m,6')
text((Hmix(5)-Hoffset_mix)/mdot_fuel/1e6,Tmix(5),'m,7')
text((Hmix(6)-Hoffset_mix)/mdot_fuel/1e6,Tmix(6),'m,8')
text((Hwater(1)-Hoffset_water)/mdot_fuel/1e6,Twater(1),'w,1')
text((Hwater(2)-Hoffset_water)/mdot_fuel/1e6,Twater(2),'w,2')
text((Hwater(3)-Hoffset_water)/mdot_fuel/1e6,Twater(3),'w,3')
text((Hwater(4)-Hoffset_water)/mdot_fuel/1e6,Twater(4),'w,4')
text((Hwater(5)-Hoffset_water)/mdot_fuel/1e6,Twater(5),'w,5')
text((Hwater(6)-Hoffset_water)/mdot_fuel/1e6,Twater(6),'w,6')
text((Hwater(7)-Hoffset_water)/mdot_fuel/1e6,Twater(7),'w,7')
text((Hwater(8)-Hoffset_water)/mdot_fuel/1e6,Twater(8),'w,8')

text(-29,1460,sprintf('Gas Turbine Power:  %.2f MW',W_net_gas_turbine*1e-6))
text(-29,1300,sprintf('Steam Turbine Power:  %.2f MW',W_net_steam_turbine*1e-6))

xlabel("Fuel-Specific Enthalpy difference (MJ/kg-fuel)");
ylabel("Temperature (K)");

hold off
plotfixer
legend('off')