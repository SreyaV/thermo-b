function dPdr = dPdr_crT(c,r,T)
% Return the derivative of P wrt r (Pa/kg/m3) for r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.

f = 1;
for i=1:1:length(c)
    f = f + c(i)*2*d*ard_idt(i,d,t) + c(i)*d*d*ardd_idt(i,d,t);
end

% Add the excess.

f = f + 2*d*aed_cdt(c,d,t) + d*d*aedd_cdt(c,d,t);

% Scale.

dPdr = R_c(c)*T*f;