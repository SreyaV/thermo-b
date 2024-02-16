function Sat = Saturation_iT(ispecies,T)
% Return the saturation pressure (Pa) and density (kg/m3) for any species i 
% at any given T (K).

rmax = 100;
rmin = 0.1;
steps = 2000;
dr = (rmax-rmin)/steps;
% Preallocate storage...
Prisotherm  = zeros(steps+1);
risotherm   = zeros(steps+1);
Chem_Potentials = zeros(steps+1);

% First generate P, r and mu arrays
i=1;
for r=rmin:dr:rmax+dr
    Prisotherm(i) = P_irT(ispecies,r,T);
    Chem_Potentials(i) = mu_irT(ispecies,r,T);
    risotherm(i) = r;
    i=i+1;
end

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

% Split into 2 curves and look for crossing
chem_pot1 = Chem_Potentials(1:i1);
chem_pot3 = Chem_Potentials(i2:end);
pressure1 = Prisotherm(1:i1);
pressure3 = Prisotherm(i2:end);

i = 1;
% Reduce the search to an interval (1st curve) i11=>i1
while chem_pot1(i)<Chem_Potentials(i2)
    i = i+1;
end
i11 = i;

i=1;
% Reduce the search to an interval (2nd curve) 1(i2)=>i22
while chem_pot3(i)<Chem_Potentials(i1)
    i = i+1;
end
i22 = i;

% FInd cross over
min_value = norm([chem_pot3(1)-chem_pot1(i11) pressure1(i11)-pressure3(1)]);
sat_P = pressure1(i11);
index = i11;
for j=i11:i1 %Looping through curve 1
    
    for l=1:i22 %Looping through curve 2
            
           distance = norm([chem_pot3(l)-chem_pot1(j) pressure3(l)-pressure1(j)]);
           if distance < min_value
               min_value = distance;
               sat_P = pressure1(j);
               index = j;
           end
    end
end
Sat.P = sat_P;
Sat.r = risotherm(index);
Sat.mu = chem_pot1(index);
end