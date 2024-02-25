function w = w_crT(c,r,T)
% Return speed of sound (m/s) for given r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

% Reduce the variables.
t = Tred_c(c)/T;
d = r/rred_c(c);

% Need the heat capacities.
Cp = Cp_crT(c,r,T);
Cv = Cv_crT(c,r,T);

% Make an ideal solution.
f = 1;
for i=1:1:length(c)
    f = f + c(i)*(...
        2*d*ard_idt(i,d,t) + d*d*ardd_idt(i,d,t)...
        );
end

% Add the excess.
f = f + 2*d*aed_cdt(c,d,t) + d*d*aedd_cdt(c,d,t);

% Scale.
w  = sqrt(R_c(c)*T*(Cp/Cv)*f);
