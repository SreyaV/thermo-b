%% Q2b table generation

clear all
% Set up the basic storage and load the FR files.
Setup_Props_i;

for ispecies = 1:6
    ispecies
    T = Ttrip_i(ispecies);
    Sat = saturation_iT(ispecies,T);
    pressure3(ispecies) = Sat.P;
    density3(ispecies)  = Sat.r;
end


%% Table
Setup_Props_i;
fprintf('\n')
fprintf('Triple Point of different Species \n')
fprintf('\n')
fprintf('Species      T(K)   Pressure (kPa)   Pressure (kPa)   Density(kg/m3)   Density(kg/m3) \n')
fprintf('                      Numerical       Experimental      Numerical       Experimental      \n')
fprintf('----------------------------------------------------------------------------------------------\n')
fprintf('N2      %10.3f %10.3f       %10.3f        %10.3f         %10.3f   \n',Ttrip_i(1),pressure3(1)/1e3, Ptrip_i(1)/1e3,density3(1),rftrip_i(1))
fprintf('O2      %10.3f %10.3f       %10.3f        %10.3f         %10.3f   \n',Ttrip_i(2),pressure3(2)/1e3, Ptrip_i(2)/1e3,density3(2),rftrip_i(2))
fprintf('Ar      %10.3f %10.3f       %10.3f        %10.3f         %10.3f   \n',Ttrip_i(3),pressure3(3)/1e3, Ptrip_i(3)/1e3,density3(3),rftrip_i(3))
fprintf('CO2     %10.3f %10.3f       %10.3f        %10.3f         %10.3f   \n',Ttrip_i(4),pressure3(4)/1e3, Ptrip_i(4)/1e3,density3(4),rftrip_i(4))
fprintf('H2O     %10.3f %10.3f       %10.3f        %10.3f         %10.3f   \n',Ttrip_i(5),pressure3(5)/1e3, Ptrip_i(5)/1e3,density3(5),rftrip_i(5))
fprintf('nH2     %10.3f %10.3f       %10.3f        %10.3f         %10.3f   \n',Ttrip_i(6),pressure3(6)/1e3, Ptrip_i(6)/1e3,density3(6),rftrip_i(6))


%%

function Sat = saturation_iT(ispecies,T)
% Return the saturation pressure (Pa) and density (kg/m3) for any species i 
% at any given T (K).
Setup_Props_i

rmax = rupper_i(ispecies);
rmin = rgtrip_i(ispecies)/2;
steps = 20000;
dr = (rmax-rmin)/steps;

% Preallocate storage...
Prisotherm  = zeros(steps+1,1);
risotherm   = zeros(steps+1,1);
Chem_Potentials = zeros(steps+1,1);

% First generate P, r and mu arrays
i=1;
for r=rmin:dr:rmax+dr
    Prisotherm(i) = P_irT(ispecies,r,T);
    Chem_Potentials(i) = mu_irT(ispecies,r,T);
    risotherm(i) = r;
    i=i+1;
end

% % Checking the curve
% scatter(Prisotherm/(1e6),Chem_Potentials/(1e6))
% plotfixer
% xlabel("P (MPa)")
% ylabel("Chem Pot (MJ/kmol)")
% legend("off")


% Cross finder

% Find turnover points i1 and i2
i = 1;
while Chem_Potentials(i+1)>Chem_Potentials(i)
    i=i+1;
end
i1 = i;

while Chem_Potentials(i+1)<Chem_Potentials(i)
    i=i+1;
end
i2 = i;

%Debugging for near triple points curves 
if Chem_Potentials(i2)<100*Chem_Potentials(i1)
    while Chem_Potentials(i+1)>=Chem_Potentials(i)
        i=i+1;
    end
    while Chem_Potentials(i+1)<=Chem_Potentials(i)
        i=i+1;
    end
    i2 = i;
end

% Split into 2 curves and look for crossing
chem_pot1 = Chem_Potentials(1:i1);
chem_pot3 = Chem_Potentials(i2:end);
pressure1 = Prisotherm(1:i1);
pressure3 = Prisotherm(i2:end);

i = 1;
% Reduce the search to an interval (1st curve) i11=>i1
while chem_pot1(i+1)<Chem_Potentials(i2)
    i = i+1;
end
i11 = i;

i=1;
% Reduce the search to an interval (2nd curve) 1(i2)=>i22
while chem_pot3(i)<Chem_Potentials(i1)
    if i+i2+1<steps
        i = i+1;
    else
        break
    end
end
i22 = i;


% Find cross over points
min_value = norm([chem_pot3(1)-chem_pot1(i11) pressure1(i11)-pressure3(1)]);
for j=i11:i1 %Looping through curve 1
    
    for l=1:i22 %Looping through curve 2
            
           distance = norm([chem_pot3(l)-chem_pot1(j) pressure3(l)-pressure1(j)]);
           if distance < min_value
               min_value = distance;
               index1 = j;
               index3 = l;
           end
    end
end


%Take the linear crossing between the 2 curves : a1x + b1 = a3x + b3
% With x = pressure, ,y = chemical potential
% This gives x = (b1-b3)/(a3-a1)
% And also y = a1x + b1
a1 = (chem_pot1(index1+1)-chem_pot1(index1))/ (pressure1(index1+1)-pressure1(index1));
a3 = (chem_pot3(index3+1)-chem_pot3(index3))/ (pressure3(index3+1)-pressure3(index3));
b1 = chem_pot1(index1)-a1*pressure1(index1);
b3 = chem_pot3(index3)-a3*pressure3(index3);

Sat.P = (b1-b3) / (a3-a1);
Sat.mu = a1 * Sat.P + b1;

Sat.r  = risotherm(index1)+risotherm(index3+i2);

end













