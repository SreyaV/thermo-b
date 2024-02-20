function [Psat rf rg] = Saturation_iT(ispecies,T)
% Return the saturation pressure (Pa) and density (kg/m3) for any species i 
% at any given T (K).
Setup_Props_i

rmax = rupper_i(ispecies);
rmin = rgtrip_i(ispecies)/2;
steps = 2000;  % CHANGE THIS NEAR TRIPLE POINT
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
% axis([1.255 1.26 -0.625 -0.624])


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

% Debugging for near triple points curves 
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

Psat = (b1-b3) / (a3-a1);
%Sat.mu = a1 * Sat.P + b1;

rg = risotherm(index1);
rf =risotherm(index3+i2);

end