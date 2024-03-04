function d2PdrdT = d2PdrdT_crT(c,r,T)
% Return the mixed derivative of P wrt r and T for r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.

f = 1;
for i=1:1:length(c)
    f = f + c(i)*(2*d*ard_idt(i,d,t) + d*d*ardd_idt(i,d,t)...
        - 2*t*d*ardt_idt(i,d,t) - t*d*d*arddt_idt(i,d,t));
end

% Add the excess.

f = f + 2*d*aedt_cdt(c,d,t) + d*d*aeddt_cdt(c,d,t);

% Scale.

d2PdrdT = R_c(c)*f;