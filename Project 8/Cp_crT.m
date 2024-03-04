function Cp = Cp_crT(c,r,T)
% Return the heat capacity (J/kg-K) for r (kg/m3) and T (K);
% C.F. Edwards, 2-16-10

t = Tred_c(c)/T;
d = r/rred_c(c);

% Need Cv first.
Cv = Cv_crT(c,r,T);

% Make an ideal solution.
fN = 1;
fD = 1;
for i=1:1:length(c)
    fN = fN + c(i)*(...
        d*ard_idt(i,d,t) - d*t*ardt_idt(i,d,t)...
        );
    fD = fD + c(i)*(...
        2*d*ard_idt(i,d,t) + d*d*ardd_idt(i,d,t)...
        );
end

% Add the excess.
fN = fN + d*aed_cdt(c,d,t) + d*t*aedt_cdt(c,d,t);
fD = fD + 2*d*aed_cdt(c,d,t) + d*d*aedd_cdt(c,d,t);

% Assemble and scale.
Cp = Cv + R_c(c)*(fN^2/fD);
