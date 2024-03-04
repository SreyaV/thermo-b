function d3Pdr2dT = d3Pdr2dT_crT(c,r,T)
% Return the third mixed derivative of P wrt r,T (Pa/kg/m3) for r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.

f = 0;
for i=1:1:length(c)
    f = f + c(i)*(2*ard_idt(i,d,t) + 2*d*ardd_idt(i,d,t)...
        + 2*d*ardd_idt(i,d,t) + d*d*arddd_idt(i,d,t)...
        - 2*t*ardt_idt(i,d,t) - 2*t*d*arddt_idt(i,d,t)...
        - 2*t*d*arddt_idt(i,d,t) - t*d*d*ardddt_idt(i,d,t));
end

% Add the excess.

f = f + 2*aedt_cdt(c,d,t) + 2*d*aeddt_cdt(c,d,t)...
    + 2*d*aeddt_cdt(c,d,t) + d*d*aedddt_cdt(c,d,t);

% Scale.

d3Pdr2dT = R_c(c)*f/rred_c(c);