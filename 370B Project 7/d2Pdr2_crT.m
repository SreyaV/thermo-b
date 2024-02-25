function d2Pdr2 = d2Pdr2_crT(c,r,T)
% Return the second derivative of P wrt r (Pa/kg/m3) for r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.

f = 0;
for i=1:1:length(c)
    f = f + c(i)*2*d*ard_idt(i,d,t) + ...
        c(i)*4*d*d*ardd_idt(i,d,t) + ...
        c(i)*d*d*d*arddd_idt(i,d,t);
end

% Add the excess.

f = f + 2*d*aed_cdt(c,d,t) + 4*d*d*aedd_cdt(c,d,t) + d*d*d*aeddd_cdt(c,d,t);

% Scale.

d2Pdr2 = R_c(c)*T*f/r;