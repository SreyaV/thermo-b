o2c=0:res*10:1;
w2c=0:res:3;

%% 1 TEMP
[x, y]=contour(o2c,w2c,temperatures);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('Temperature')
plotfixer
legend('off')

%% 2 EXERGY
exergy_eff=exergy_eff*100;
[x, y]=contour(o2c,w2c,exergy_eff,[98, 96, 94, 92, 90, 88]);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('Exergy Efficiency(%)')
plotfixer
legend('off')

%% 3 Cold Gas Eff
coldGas_eff=100*coldGas_eff;
[x, y]=contour(o2c,w2c,coldGas_eff);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('Cold Gas Efficiency(%)')
plotfixer
legend('off')

%% 4 H2 Yield
[x, y]=contour(o2c,w2c,H2Yield);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('H2 Mole Fractions')
plotfixer
legend('off')

%% 5 CO Yield
[x, y]=contour(o2c,w2c,COYield);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('CO Mole Fractions')
plotfixer
legend('off')

%% 6 CO2 Yield
[x, y]=contour(o2c,w2c,CO2Yield);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('CO2 Mole Fractions')
plotfixer
legend('off')

%% 7 O2 Yield
[x, y]=contour(o2c,w2c,H2OYield);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('O2 Mole Fractions')
plotfixer
legend('off')

%% 8 CH4 Yield
[x, y]=contour(o2c,w2c,CH4Yield,[1e-6, 1e-5, 1e-4 ,1e-3, 1e-2, 0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('CH4 Mole Fractions')
plotfixer
legend('off')

%% 9 N2 Yield
[x, y]=contour(o2c,w2c,N2Yield);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
hold off
title('N2 Mole Fractions')
plotfixer
legend('off')

%% 10 H2O Yield
[x, y]=contour(o2c,w2c,H2OYield);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
hold off
title('H2O Mole Fractions')
plotfixer
legend('off')

%% 11 Syngas Yield
[x, y]=contour(o2c,w2c,syngasYield);
clabel(x,y)
xlabel('Oxygen/Carbon Molar Feed Ratio')
ylabel('Water/Carbon Molar Feed Ratio')
title('Syngas Yield')
plotfixer
legend('off')