%% Q2b table generation

% Set up the basic storage and load the FR files.
Setup_Props_i;

for ispecies = 1:6
    ispecies
    T = Ttrip_i(ispecies);
    [Psat rf rg] = Saturation_iT(ispecies,T);
    pressure3(ispecies) = Psat;
    density3f(ispecies)  = rf;
    density3g(ispecies)  = rg;
end


%% Table
Setup_Props_i;
fprintf('\n')
fprintf('Triple Point of different Species (Resolution in density: 100k points) \n')
fprintf('\n')
fprintf('Species      T(K)   Pressure (kPa)   Pressure (kPa)   Density Liq(kg/m3)   Density Liq(kg/m3)   Density Vap(kg/m3)   Density Vap(kg/m3)  \n')
fprintf('                      Numerical       Experimental       Numerical            Experimental         Numerical            Experimental     \n')
fprintf('-----------------------------------------------------------------------------------------------------------------------------------------\n')
fprintf('N2      %10.3f  %10.6f        %10.6f         %10.6f            %10.6f          %10.6f         %10.6f\n',Ttrip_i(1),pressure3(1)/1e3, Ptrip_i(1)/1e3,density3f(1),rftrip_i(1),density3g(1),rgtrip_i(1))
fprintf('O2      %10.3f  %10.6f        %10.6f        %10.6f           %10.6f          %10.6f         %10.6f\n',Ttrip_i(2),pressure3(2)/1e3, Ptrip_i(2)/1e3,density3f(2),rftrip_i(2),density3g(2),rgtrip_i(2))
fprintf('Ar      %10.3f  %10.6f        %10.6f        %10.6f           %10.6f          %10.6f         %10.6f\n',Ttrip_i(3),pressure3(3)/1e3, Ptrip_i(3)/1e3,density3f(3),rftrip_i(3),density3g(3),rgtrip_i(3))
fprintf('CO2     %10.3f  %10.6f        %10.6f        %10.6f           %10.6f          %10.6f         %10.6f\n',Ttrip_i(4),pressure3(4)/1e3, Ptrip_i(4)/1e3,density3f(4),rftrip_i(4),density3g(4),rgtrip_i(4))
fprintf('H2O     %10.3f  %10.6f        %10.6f         %10.6f            %10.6f          %10.6f         %10.6f\n',Ttrip_i(5),pressure3(5)/1e3, Ptrip_i(5)/1e3,density3f(5),rftrip_i(5),density3g(5),rgtrip_i(5))
fprintf('nH2     %10.3f  %10.6f        %10.6f         %10.6f            %10.6f          %10.6f         %10.6f\n',Ttrip_i(6),pressure3(6)/1e3, Ptrip_i(6)/1e3,density3f(6),rftrip_i(6),density3g(6),rgtrip_i(6))


