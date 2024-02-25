function d3Pdr3 = d3Pdr3_crT(c,r,T)
% Return the third derivative of P wrt r (Pa/kg/m3) for r (kg/m3) and T (K).
% C.F. Edwards, 2-16-10

t = Tred_c(c)/T;
d = r/rred_c(c);

% Make an ideal solution.

f = 0;
for i=1:1:length(c)
    f = f + c(i)*(6*d*d*ardd_idt(i,d,t) + ...
        6*d*d*d*arddd_idt(i,d,t) + ...
        d*d*d*d*ardddd_idt(i,d,t));
end

% Add the excess.

f = f + 6*d*d*aedd_cdt(c,d,t)...
    + 6*d*d*d*aeddd_cdt(c,d,t)...
    + d*d*d*d*aedddd_cdt(c,d,t);

% Scale.

d3Pdr3 = R_c(c)*T*f/r/r;