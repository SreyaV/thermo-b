function d3Pdr3 = d3Pdr3_rT(i,r,T)
% Return the third derivative of P wrt r (Pa/kg/m3) for r (kg/m3) and T (K).
% C.F. Edwards, 2-3-10

global Tcrit_i rcrit_i R_i

t = Tcrit_i(i)/T;
d = r/rcrit_i(i);

% Make an ideal solution.

f = 6*d*d*ardd_idt(i,d,t) + 6*d*d*d*arddd_idt(i,d,t) + d*d*d*d*ardddd_idt(i,d,t);

% Scale.

d3Pdr3 = R_i(i)*T*f/r/r;
